"""Loci tab: the per-locus table, explicit subset analysis, estimates.

State distinctions (docs/dev/SCIENTIFIC_CONTRACT_locus_features.md):
browsing/sorting/selecting here changes NOTHING analytically; only the
explicit "Analyze selected loci" / "Compare with eligible" actions
compute a (clearly labeled) subset aggregate — from the original
per-copy values, never by re-estimation. The eligible population's
aggregates are never touched.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from PySide6.QtCore import QAbstractTableModel, QSortFilterProxyModel, Qt
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSplitter,
    QTableView,
    QVBoxLayout,
    QWidget,
)

from ...core.loci import export_selection, import_selection, subset_result
from ...core.loci.estimates import EstimateError, load_estimates, match_estimates
from ...core.loci.selection import SelectionError
from ...core.loci.subset import SubsetError

_ID_COL = "locus_id"

#: identity/class columns shown before per-signal columns
_FIXED = [
    ("locus_id", "Locus ID"),
    ("completeness", "Completeness"),
    ("consensus_coverage", "Cons. coverage"),
    ("genome_length", "Length (bp)"),
    ("rm_id", "RM ID"),
    ("rm_id_fragments", "RM frags"),
]


class LocusTableModel(QAbstractTableModel):
    """Read-only model over the locus table DataFrame."""

    def __init__(self, table: pd.DataFrame, est_meta: dict | None = None):
        super().__init__()
        self._df = table.reset_index(drop=False).rename(
            columns={"index": "copy_idx"}
        )
        cols = [(k, h) for k, h in _FIXED if k in self._df.columns]
        for c in table.columns:
            if c.startswith("mean_"):
                cols.append((c, f"mean {c[5:]}"))
            elif c.startswith("obs_bins_"):
                cols.append((c, f"obs bins {c[9:]}"))
            elif c.startswith("peaks_"):
                cols.append((c, f"peaks {c[6:]}"))
            elif c.startswith("est_"):
                unit = ""
                if est_meta and c == "est_total":
                    unit = f" ({est_meta.get('units', '?')})"
                cols.append((c, c.replace("est_", "estimate ") + unit))
        self._cols = cols
        self._tooltips = {
            "locus_id": "Stable locus ID: chrom:start-end(strand), "
                        "0-based half-open (machine key used by selection "
                        "and estimate files). Cursor readouts elsewhere "
                        "are 1-based.",
            "rm_id_fragments": "Analyzed records sharing this RepeatMasker "
                               "ID (fragment linkage; records are never "
                               "merged).",
        }

    # Qt model API ---------------------------------------------------------
    def rowCount(self, parent=None):  # noqa: N802
        return 0 if parent is not None and parent.isValid() else len(self._df)

    def columnCount(self, parent=None):  # noqa: N802
        return 0 if parent is not None and parent.isValid() else len(self._cols)

    def headerData(self, section, orientation, role=Qt.DisplayRole):  # noqa: N802
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self._cols[section][1]
        if orientation == Qt.Horizontal and role == Qt.ToolTipRole:
            return self._tooltips.get(self._cols[section][0], "")
        return None

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        key = self._cols[index.column()][0]
        val = self._df.iloc[index.row()][key]
        if role == Qt.UserRole:  # raw value for numeric sorting
            if isinstance(val, float) and np.isnan(val):
                return float("-inf")  # sort missing last/first consistently
            return val
        if role == Qt.DisplayRole:
            if isinstance(val, float):
                if np.isnan(val):
                    # missing evidence — explicitly not a zero
                    return "unavailable"
                if key == "consensus_coverage":
                    return f"{val:.1%}"
                return f"{val:.4g}"
            return str(val)
        if role == Qt.ToolTipRole and isinstance(val, float) and np.isnan(val):
            return ("No evidence for this locus in this column "
                    "(missing/unsupported — NOT an observed zero).")
        return None

    # helpers --------------------------------------------------------------
    def locus_id_at(self, source_row: int) -> str:
        return str(self._df.iloc[source_row][_ID_COL])

    def copy_position_at(self, source_row: int) -> int:
        """Positional row in result.copies order."""
        return int(source_row)

    def row_of_locus(self, locus_id: str) -> int | None:
        hits = np.where(self._df[_ID_COL].to_numpy() == locus_id)[0]
        return int(hits[0]) if hits.size else None


class LociTab(QWidget):
    """The Loci panel + explicit subset actions + comparison canvas."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.result = None
        self.model: LocusTableModel | None = None
        self.subset = None  # applied SubsetResult | None
        self._estimates_file = ""
        self._inspector = None
        # loaded-project state carried until the next run produces a result
        self._pending_selected: list[str] = []
        self._pending_applied = False

        lay = QVBoxLayout(self)

        top = QHBoxLayout()
        self.search_edit = QLineEdit()
        self.search_edit.setPlaceholderText(
            "Filter loci (ID, chromosome, completeness…)")
        self.search_edit.setClearButtonEnabled(True)
        top.addWidget(self.search_edit, 1)
        self.inspect_btn = QPushButton("Inspect…")
        self.inspect_btn.setToolTip(
            "Open the Locus Inspector for the focused locus "
            "(also: double-click a row)")
        top.addWidget(self.inspect_btn)
        self.estimates_btn = QPushButton("Attach estimates…")
        self.estimates_btn.setToolTip(
            "Attach an upstream per-locus estimate table "
            "(tescope-locus-estimates TSV). Locus-level table context "
            "only — never drawn as positional signal.")
        top.addWidget(self.estimates_btn)
        lay.addLayout(top)

        self.table = QTableView()
        self.table.setSelectionBehavior(QTableView.SelectRows)
        self.table.setSelectionMode(QTableView.ExtendedSelection)
        self.table.setSortingEnabled(True)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)

        self.proxy = QSortFilterProxyModel(self)
        self.proxy.setSortRole(Qt.UserRole)
        self.proxy.setFilterCaseSensitivity(Qt.CaseInsensitive)
        self.proxy.setFilterKeyColumn(-1)  # all columns
        self.table.setModel(self.proxy)

        actions = QHBoxLayout()
        self.analyze_btn = QPushButton("Analyze selected loci")
        self.analyze_btn.setToolTip(
            "EXPLICIT action: aggregate only the selected loci "
            "(recomputed from the original per-copy values; the "
            "eligible population's profile is untouched).")
        self.compare_btn = QPushButton("Compare with eligible")
        self.compare_btn.setToolTip(
            "EXPLICIT action: draw the selected-locus aggregate on top "
            "of the eligible population's aggregate, with both "
            "denominators shown.")
        self.clear_btn = QPushButton("Clear subset")
        self.export_sel_btn = QPushButton("Export selection…")
        self.import_sel_btn = QPushButton("Import selection…")
        self.export_fig_btn = QPushButton("Export figure…")
        for b in (self.analyze_btn, self.compare_btn, self.clear_btn,
                  self.export_sel_btn, self.import_sel_btn,
                  self.export_fig_btn):
            actions.addWidget(b)
        actions.addStretch(1)
        lay.addLayout(actions)

        self.status = QLabel("Run an analysis to list its loci.")
        self.status.setObjectName("hintLabel")
        self.status.setWordWrap(True)

        split = QSplitter(Qt.Vertical)
        split.addWidget(self.table)
        from .results import _CanvasHolder  # shared canvas container

        self.canvas = _CanvasHolder(result_provider=lambda: self.result)
        split.addWidget(self.canvas)
        split.setSizes([380, 240])
        lay.addWidget(split, 1)
        lay.addWidget(self.status)

        # wiring
        self.search_edit.textChanged.connect(self.proxy.setFilterFixedString)
        self.analyze_btn.clicked.connect(lambda: self._apply_subset(False))
        self.compare_btn.clicked.connect(lambda: self._apply_subset(True))
        self.clear_btn.clicked.connect(self._clear_subset)
        self.export_sel_btn.clicked.connect(self._export_selection)
        self.import_sel_btn.clicked.connect(self._import_selection)
        self.estimates_btn.clicked.connect(self._attach_estimates)
        self.export_fig_btn.clicked.connect(self._export_figure)
        self.inspect_btn.clicked.connect(self._open_inspector)
        self.table.doubleClicked.connect(lambda _i: self._open_inspector())
        self._set_enabled(False)

    # -- result lifecycle ---------------------------------------------------

    def set_result(self, result) -> None:
        self.result = result
        self.subset = None
        if self._inspector is not None:
            self._inspector.close()
            self._inspector = None
        if result is None:
            self.model = None
            self.proxy.setSourceModel(None)
            self._set_enabled(False)
            self.status.setText("Run an analysis to list its loci.")
            return
        if result.config.loci.estimates_file:
            self._estimates_file = result.config.loci.estimates_file
        est_meta = (result.locus_estimates.to_dict()
                    if result.locus_estimates is not None else None)
        self.model = LocusTableModel(result.locus_table(), est_meta)
        self.proxy.setSourceModel(self.model)
        self._set_enabled(True)
        self.canvas.set_figure(None)
        self.table.selectionModel().selectionChanged.connect(
            lambda *_: self._update_status())
        # restore a persisted selection (stale IDs are reported)
        persisted = list(result.config.loci.selected or [])
        if persisted:
            stale = self._select_ids(persisted)
            if stale:
                QMessageBox.information(
                    self, "Selection partially restored",
                    f"{len(stale)} saved locus ID(s) are not in the "
                    "current eligible population (filters or annotation "
                    "changed):\n" + "\n".join(stale[:10]))
            if result.config.loci.subset_applied:
                self._apply_subset(compare=True)
        self._update_status()

    def _set_enabled(self, on: bool) -> None:
        for w in (self.search_edit, self.table, self.analyze_btn,
                  self.compare_btn, self.clear_btn, self.export_sel_btn,
                  self.import_sel_btn, self.estimates_btn,
                  self.export_fig_btn, self.inspect_btn):
            w.setEnabled(on)

    # -- selection helpers --------------------------------------------------

    def selected_locus_ids(self) -> list[str]:
        if self.model is None:
            return []
        rows = {self.proxy.mapToSource(i).row()
                for i in self.table.selectionModel().selectedRows()}
        return [self.model.locus_id_at(r) for r in sorted(rows)]

    def _select_ids(self, ids: list[str]) -> list[str]:
        """Select rows for the given IDs; returns the stale ones."""
        sel = self.table.selectionModel()
        sel.clearSelection()
        stale = []
        from PySide6.QtCore import QItemSelection, QItemSelectionModel

        selection = QItemSelection()
        for lid in ids:
            row = self.model.row_of_locus(lid)
            if row is None:
                stale.append(lid)
                continue
            src = self.model.index(row, 0)
            pi = self.proxy.mapFromSource(src)
            if pi.isValid():
                selection.select(pi, self.proxy.index(
                    pi.row(), self.proxy.columnCount() - 1))
        sel.select(selection, QItemSelectionModel.Select | QItemSelectionModel.Rows)
        return stale

    def _update_status(self) -> None:
        if self.result is None:
            return
        n = len(self.result.copies)
        k = len(self.selected_locus_ids())
        parts = [f"Eligible population: {n:,} loci", f"selected: {k:,}"]
        if self.subset is not None:
            parts.append(f"APPLIED {self.subset.label}")
            if self.subset.invalid_ids:
                parts.append(
                    f"{len(self.subset.invalid_ids)} saved ID(s) invalid")
        else:
            parts.append("no subset applied — selection changes nothing "
                         "until applied")
        if self.result.locus_estimates is not None:
            est = self.result.locus_estimates
            rep = self.result.estimate_report
            parts.append(
                f"estimates: {est.assay}, {est.units}"
                + (f", {rep.summary()}" if rep is not None else ""))
        self.status.setText("  |  ".join(parts))

    # -- explicit subset actions --------------------------------------------

    def _apply_subset(self, compare: bool) -> None:
        ids = self.selected_locus_ids()
        if self.result is None:
            return
        try:
            self.subset = subset_result(self.result, ids)
        except SubsetError as e:
            QMessageBox.warning(self, "Cannot analyze subset", str(e))
            return
        # persist the explicit state for reproducibility (project save)
        self.result.config.loci.selected = list(self.subset.locus_ids)
        self.result.config.loci.subset_applied = True
        from ...plotting.profiles import build_subset_comparison
        from ...plotting.style import publication_style

        with publication_style():
            fig = build_subset_comparison(self.result, self.subset,
                                          compare=compare)
        self.canvas.set_figure(fig)
        self._update_status()

    def _clear_subset(self) -> None:
        self.subset = None
        if self.result is not None:
            self.result.config.loci.subset_applied = False
        self.canvas.set_figure(None)
        self._update_status()

    # -- selection files ----------------------------------------------------

    def _export_selection(self) -> None:
        ids = self.selected_locus_ids()
        if not ids:
            QMessageBox.information(self, "Nothing selected",
                                    "Select loci to export first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export locus selection", "selection.yaml",
            "TEscope locus selection (*.yaml *.yml)")
        if not path:
            return
        try:
            export_selection(self.result, ids, path)
        except SelectionError as e:
            QMessageBox.warning(self, "Export failed", str(e))

    def _import_selection(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Import locus selection", "",
            "TEscope locus selection (*.yaml *.yml)")
        if not path:
            return
        try:
            valid, invalid, warn = import_selection(path, self.result)
        except (SelectionError, OSError) as e:
            QMessageBox.warning(self, "Import failed", str(e))
            return
        self._select_ids(valid)
        msgs = [f"{len(valid)} locus/loci selected."]
        if warn:
            msgs.append(warn)
        if invalid:
            msgs.append(
                f"{len(invalid)} ID(s) not in the current eligible "
                "population:\n" + "\n".join(invalid[:10]))
        if warn or invalid:
            QMessageBox.information(self, "Selection imported",
                                    "\n".join(msgs))
        self._update_status()

    # -- estimates ----------------------------------------------------------

    def _attach_estimates(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Attach upstream locus estimates", "",
            "tescope-locus-estimates (*.tsv *.txt)")
        if not path:
            return
        try:
            est = load_estimates(path)
            table, report = match_estimates(
                est, self.result.locus_table().drop(
                    columns=[c for c in self.result.locus_table().columns
                             if c.startswith("est_")]),
                self.result.config.assembly)
        except (EstimateError, OSError) as e:
            QMessageBox.warning(self, "Estimate import failed", str(e))
            return
        self.result.locus_estimates = est
        self.result.estimate_report = report
        self.result._locus_table = table
        self.result.config.loci.estimates_file = path
        self._estimates_file = path
        selected = self.selected_locus_ids()
        self.model = LocusTableModel(table, est.to_dict())
        self.proxy.setSourceModel(self.model)
        self.table.selectionModel().selectionChanged.connect(
            lambda *_: self._update_status())
        self._select_ids(selected)
        summary = report.summary()
        QMessageBox.information(
            self, "Estimates attached",
            f"{est.assay} estimates ({est.units}"
            + (f", tool: {est.tool}" if est.tool else "") + ").\n"
            f"{summary}\n\n"
            "Table context only: estimates are never drawn as "
            "positional signal, and selecting loci never re-fits or "
            "reassigns them.")
        self._update_status()

    # -- figure export ------------------------------------------------------

    def _export_figure(self) -> None:
        if self.canvas.canvas is None or self.subset is None:
            QMessageBox.information(
                self, "No subset figure",
                "Apply a subset first (Analyze selected loci / "
                "Compare with eligible).")
            return
        from .results import _EXPORT_FORMATS, _resolve_export

        path, chosen = QFileDialog.getSaveFileName(
            self, "Export subset figure", "tescope_subset.png",
            ";;".join(f"{f.upper()} (*.{f})" for f in _EXPORT_FORMATS))
        if not path:
            return
        stem, fmt = _resolve_export(path, chosen)
        from ...export.figures import save_figure

        try:
            written = save_figure(self.canvas.canvas.figure, stem, [fmt],
                                  result=self.result)
        except Exception as e:  # noqa: BLE001 - surface any export failure
            QMessageBox.warning(self, "Export failed", str(e))
            return
        self.status.setText(
            f"Figure written to {written[0]} (+ provenance sidecar)")

    # -- inspector ----------------------------------------------------------

    def _proxy_order_ids(self) -> list[str]:
        return [
            self.model.locus_id_at(self.proxy.mapToSource(
                self.proxy.index(r, 0)).row())
            for r in range(self.proxy.rowCount())
        ]

    def _open_inspector(self) -> None:
        if self.result is None or self.model is None:
            return
        idx = self.table.selectionModel().currentIndex()
        ids = self._proxy_order_ids()
        if not ids:
            return
        current = ids[idx.row()] if idx.isValid() and idx.row() < len(ids) \
            else ids[0]
        from .inspector import LocusInspector

        if self._inspector is None:
            self._inspector = LocusInspector(self)
            self._inspector.highlight_requested.connect(
                self._on_highlight_requested)
        self._inspector.show_locus(self.result, ids, current)
        self._inspector.show()
        self._inspector.raise_()

    def _on_highlight_requested(self, locus_id: str) -> None:
        # delegate to ResultsView (parent chain) to switch to the heatmap
        parent = self.parent()
        while parent is not None and not hasattr(
                parent, "highlight_locus_in_heatmap"):
            parent = parent.parent()
        if parent is not None:
            parent.highlight_locus_in_heatmap(locus_id)

    # -- config round-trip --------------------------------------------------

    def apply_loci_config(self, lcfg) -> None:
        """Project-load path (no result yet): stash for the next run."""
        self._estimates_file = lcfg.estimates_file
        self._pending_selected = list(lcfg.selected or [])
        self._pending_applied = bool(lcfg.subset_applied)

    def collect_into(self, cfg) -> None:
        cfg.loci.estimates_file = self._estimates_file
        if self.result is not None and self.model is not None:
            cfg.loci.selected = (list(self.subset.locus_ids) if self.subset
                                 else self.selected_locus_ids())
            cfg.loci.subset_applied = self.subset is not None
        else:
            cfg.loci.selected = list(self._pending_selected)
            cfg.loci.subset_applied = self._pending_applied
