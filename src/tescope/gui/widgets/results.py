"""Results area: profile figure, heatmap, QC table, export actions."""

from __future__ import annotations

from pathlib import Path

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from ...core.analysis import AnalysisResult
from ...core.qc import qc_table
from ...core.sequence.capabilities import sequence_capabilities
from ...core.sequence.coords import RangeError, visible_base_range
from ...export.figures import save_figure
from ...export.tables import export_tables
from ...plotting.profiles import build_figure, build_heatmap
from ...plotting.style import publication_style
from .browser import ConsensusNavigator


class _ProvenanceToolbar(NavigationToolbar2QT):
    """Navigation toolbar whose Save writes TEscope's provenance sidecar.

    This closes the historically documented bypass: toolbar Save used to
    write a bare image with no ``.figure.json``. When a result provider
    is available, Save now goes through the same export path as the
    "Export figure…" button; without one it falls back to matplotlib's
    plain save.
    """

    def __init__(self, canvas, parent, result_provider=None):
        super().__init__(canvas, parent)
        self._result_provider = result_provider

    def save_figure(self, *args):  # noqa: D102 - matplotlib signature
        result = self._result_provider() if self._result_provider else None
        if result is None:
            return super().save_figure(*args)
        path, chosen = QFileDialog.getSaveFileName(
            self.parent(), "Save figure (with provenance sidecar)",
            "tescope_figure.png",
            ";;".join(f"{f.upper()} (*.{f})" for f in _EXPORT_FORMATS))
        if not path:
            return None
        stem, fmt = _resolve_export(path, chosen)
        try:
            save_figure(self.canvas.figure, stem, [fmt], result=result)
        except Exception as e:  # noqa: BLE001 - surface any export failure
            QMessageBox.warning(self.parent(), "Save failed", str(e))
        return None


class _CanvasHolder(QWidget):
    """A widget that swaps in freshly built matplotlib figures.

    Styled as a card (``canvasFrame``); the publication figure inside is
    deliberately untouched by the application theme.
    """

    def __init__(self, parent=None, result_provider=None):
        super().__init__(parent)
        self.setObjectName("canvasFrame")
        self.setAttribute(Qt.WA_StyledBackground, True)
        self._lay = QVBoxLayout(self)
        self._lay.setContentsMargins(6, 4, 6, 6)
        self.canvas: FigureCanvasQTAgg | None = None
        self.toolbar: NavigationToolbar2QT | None = None
        self._result_provider = result_provider
        self._placeholder = QLabel(
            "Select a genome and TE, add signal tracks, then Run analysis."
        )
        self._placeholder.setObjectName("hintLabel")
        self._placeholder.setAlignment(Qt.AlignCenter)
        self._lay.addWidget(self._placeholder, 1)

    def set_figure(self, fig) -> None:
        if self._placeholder is not None:
            self._lay.removeWidget(self._placeholder)
            self._placeholder.deleteLater()
            self._placeholder = None
        if self.canvas is not None:
            self._lay.removeWidget(self.toolbar)
            self._lay.removeWidget(self.canvas)
            self.toolbar.deleteLater()
            self.canvas.deleteLater()
            self.canvas = None
            self.toolbar = None
        if fig is None:
            return
        self.canvas = FigureCanvasQTAgg(fig)
        self.toolbar = _ProvenanceToolbar(
            self.canvas, self, result_provider=self._result_provider)
        self._lay.addWidget(self.toolbar)
        self._lay.addWidget(self.canvas, 1)
        self.canvas.draw_idle()



_EXPORT_FORMATS = ("png", "svg", "pdf")


def _resolve_export(path: str, chosen: str) -> tuple:
    """(stem_path, fmt) for a figure save dialog result: a recognized
    typed extension wins; otherwise the selected filter decides; dotted
    filenames keep their full stem."""
    p = Path(path)
    ext = p.suffix.lstrip(".").lower()
    if ext in _EXPORT_FORMATS:
        return p.with_name(p.name[: -len(p.suffix)]), ext
    fmt = next((f for f in _EXPORT_FORMATS if f.upper() in (chosen or "")),
               "png")
    return p, fmt


class ResultsView(QTabWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.result: AnalysisResult | None = None
        self._current_fig = None

        # Profile tab
        profile_tab = QWidget()
        lay = QVBoxLayout(profile_tab)
        lay.setContentsMargins(10, 8, 10, 10)
        lay.setSpacing(8)
        controls = QHBoxLayout()
        controls.setSpacing(8)
        controls.addWidget(QLabel("View:"))
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["stacked", "overlay", "fragments"])
        controls.addWidget(self.mode_combo)
        controls.addWidget(QLabel("Y axis:"))
        self.yaxis_combo = QComboBox()
        self.yaxis_combo.addItems(["independent", "shared"])
        controls.addWidget(self.yaxis_combo)
        controls.addWidget(QLabel("Normalize:"))
        self.norm_combo = QComboBox()
        self.norm_combo.addItems(["none", "zscore", "minmax"])
        controls.addWidget(self.norm_combo)
        self.sem_check = QCheckBox("SEM band")
        self.sem_check.setChecked(True)
        controls.addWidget(self.sem_check)
        self.reps_check = QCheckBox("Replicates")
        controls.addWidget(self.reps_check)
        controls.addStretch(1)
        style_btn = QPushButton("Style…")
        style_btn.clicked.connect(self._open_style_dialog)
        controls.addWidget(style_btn)
        export_fig = QPushButton("Export figure…")
        export_fig.clicked.connect(self._export_figure)
        controls.addWidget(export_fig)
        export_data = QPushButton("Export data…")
        export_data.clicked.connect(self._export_data)
        controls.addWidget(export_data)
        lay.addLayout(controls)

        # -- consensus browser row: navigation + sequence functions -------
        # (display-only interaction; capability-gated per current result)
        from PySide6.QtWidgets import QLineEdit

        nav_row = QHBoxLayout()
        nav_row.setSpacing(8)
        nav_row.addWidget(QLabel("Go to:"))
        self.goto_edit = QLineEdit()
        self.goto_edit.setPlaceholderText("e.g. 4127 or 4000-4500 (1-based)")
        self.goto_edit.setMaximumWidth(180)
        self.goto_edit.returnPressed.connect(self._go_to)
        nav_row.addWidget(self.goto_edit)
        self.full_btn = QPushButton("Full view")
        self.full_btn.setToolTip("Show the whole consensus (also: toolbar Home)")
        self.full_btn.clicked.connect(lambda: self.navigator.full_view())
        nav_row.addWidget(self.full_btn)
        self.seq_check = QCheckBox("Show sequence")
        self.seq_check.toggled.connect(self._toggle_sequence)
        nav_row.addWidget(self.seq_check)
        self.copy_seq_btn = QPushButton("Copy sequence")
        self.copy_seq_btn.clicked.connect(self._copy_sequence)
        nav_row.addWidget(self.copy_seq_btn)
        self.fasta_btn = QPushButton("Export FASTA…")
        self.fasta_btn.clicked.connect(self._export_fasta)
        nav_row.addWidget(self.fasta_btn)
        nav_row.addWidget(QLabel("Search:"))
        self.search_edit = QLineEdit()
        self.search_edit.setPlaceholderText("ACGTN… (N = any)")
        self.search_edit.setMaximumWidth(140)
        self.search_edit.returnPressed.connect(self._search_sequence)
        nav_row.addWidget(self.search_edit)
        self.hits_combo = QComboBox()
        self.hits_combo.setMinimumWidth(130)
        self.hits_combo.activated.connect(self._go_to_hit)
        nav_row.addWidget(self.hits_combo)
        nav_row.addStretch(1)
        lay.addLayout(nav_row)

        self.profile_canvas = _CanvasHolder(
            result_provider=lambda: self.result)
        lay.addWidget(self.profile_canvas, 1)

        # cursor readout: plain label, updated on hover (no canvas redraw)
        self.cursor_label = QLabel(
            "Hover for position/base/signal readout — wheel zooms, "
            "toolbar pans; click marks a position.")
        self.cursor_label.setObjectName("hintLabel")
        lay.addWidget(self.cursor_label)
        self.addTab(profile_tab, "Profile")

        self.navigator = ConsensusNavigator(self)
        self.navigator.cursor_info.connect(self.cursor_label.setText)
        self._search_hits = []

        for combo in (self.mode_combo, self.yaxis_combo, self.norm_combo):
            combo.currentTextChanged.connect(lambda *_: self.redraw())
        self.sem_check.toggled.connect(lambda *_: self.redraw())
        self.reps_check.toggled.connect(lambda *_: self.redraw())

        # Heatmap tab
        heat_tab = QWidget()
        hlay = QVBoxLayout(heat_tab)
        hlay.setContentsMargins(10, 8, 10, 10)
        hlay.setSpacing(8)
        hrow = QHBoxLayout()
        hrow.setSpacing(8)
        hrow.addWidget(QLabel("Track:"))
        self.heat_combo = QComboBox()
        self.heat_combo.currentTextChanged.connect(lambda *_: self._redraw_heatmap())
        hrow.addWidget(self.heat_combo, 1)
        heat_export = QPushButton("Export heatmap…")
        heat_export.clicked.connect(self._export_heatmap)
        hrow.addWidget(heat_export)
        hlay.addLayout(hrow)
        self.heat_canvas = _CanvasHolder(
            result_provider=lambda: self.result)
        hlay.addWidget(self.heat_canvas, 1)
        self.addTab(heat_tab, "Heatmap")
        self._heat_highlight: int | None = None

        # Loci tab (per-locus table + explicit subset analysis)
        from .loci import LociTab

        self.loci_tab = LociTab(self)
        self.addTab(self.loci_tab, "Loci")

        # QC tab
        self.qc_widget = QTableWidget(0, 3)
        self.qc_widget.setHorizontalHeaderLabels(["Section", "Metric", "Value"])
        self.qc_widget.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.qc_widget.verticalHeader().setVisible(False)
        self.qc_widget.setShowGrid(False)
        self.qc_widget.setAlternatingRowColors(True)
        self.addTab(self.qc_widget, "QC")

    # -- population --------------------------------------------------------

    def show_result(self, result: AnalysisResult) -> None:
        self.result = result
        # sync controls from config without triggering redraw storms
        for combo, val in (
            (self.yaxis_combo, result.config.plot.y_axis),
            (self.norm_combo, result.config.plot.normalization),
        ):
            combo.blockSignals(True)
            combo.setCurrentText(val)
            combo.blockSignals(False)
        self.heat_combo.blockSignals(True)
        self.heat_combo.clear()
        self.heat_combo.addItems(
            [t.name for t in result.continuous_tracks() if t.matrix is not None]
        )
        self.heat_combo.blockSignals(False)
        self._update_sequence_gating(result)
        self._heat_highlight = None
        self.loci_tab.set_result(result)
        self.redraw()
        self._redraw_heatmap()
        self._fill_qc(result)

    def redraw(self) -> None:
        if self.result is None:
            return
        cfg = self.result.config
        cfg.plot.y_axis = self.yaxis_combo.currentText()
        cfg.plot.normalization = self.norm_combo.currentText()
        cfg.plot.show_sem = self.sem_check.isChecked()
        cfg.plot.show_replicates = self.reps_check.isChecked()
        with publication_style():
            fig = build_figure(self.result, mode=self.mode_combo.currentText())
        self._current_fig = fig
        self.profile_canvas.set_figure(fig)
        # (re)bind interactive navigation to the fresh canvas
        self.navigator.attach(self.profile_canvas.canvas, self.result)

    def _redraw_heatmap(self) -> None:
        if self.result is None or not self.heat_combo.currentText():
            return
        name = self.heat_combo.currentText()
        track = next((t for t in self.result.continuous_tracks() if t.name == name), None)
        if track is None:
            return
        with publication_style():
            fig = build_heatmap(self.result, track,
                                highlight_row=self._heat_highlight)
        self.heat_canvas.set_figure(fig)

    def highlight_locus_in_heatmap(self, locus_id: str) -> None:
        """Display-only: mark one locus's row in the per-copy heatmap."""
        if self.result is None:
            return
        import numpy as np

        from ...core.loci.table import locus_ids

        ids = locus_ids(self.result.copies).to_numpy()
        pos = np.where(ids == locus_id)[0]
        if not pos.size:
            return
        self._heat_highlight = int(pos[0])
        self._redraw_heatmap()
        heat_index = next(
            i for i in range(self.count())
            if self.widget(i) is not None and self.tabText(i) == "Heatmap")
        self.setCurrentIndex(heat_index)

    def _fill_qc(self, result: AnalysisResult) -> None:
        qc = qc_table(result)
        self.qc_widget.setRowCount(len(qc))
        for i, r in enumerate(qc.itertuples()):
            for j, val in enumerate((r.section, r.metric, r.value)):
                item = QTableWidgetItem(str(val))
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.qc_widget.setItem(i, j, item)
        self.qc_widget.resizeColumnsToContents()

    # -- consensus browser / sequence functions ----------------------------

    def _seq_caps(self):
        res = None if self.result is None else getattr(self.result, "sequence", None)
        return sequence_capabilities(res)

    def _update_sequence_gating(self, result: AnalysisResult) -> None:
        """Capability-driven enable/disable with explanatory tooltips.

        Conditions come from the centralized sequence-capability model —
        controls are visible but disabled when unavailable, and the
        tooltip says why (same philosophy as the track Mapping controls).
        """
        caps = self._seq_caps()
        for widget, cap, tip_ok in (
            (self.seq_check, "base_display",
             "Show the aligned consensus sequence track (verified model)"),
            (self.copy_seq_btn, "copy_sequence",
             "Copy the visible consensus sequence (FASTA, 1-based range)"),
            (self.fasta_btn, "export_fasta",
             "Export the consensus sequence or visible region as FASTA"),
            (self.search_edit, "search",
             "Literal sequence search (case-insensitive, N = any base)"),
            (self.hits_combo, "search", "Search hits — select to zoom"),
        ):
            ok = getattr(caps, cap)
            widget.setEnabled(ok)
            widget.setToolTip(tip_ok if ok else caps.why_not(cap))
        self.seq_check.blockSignals(True)
        self.seq_check.setChecked(
            caps.verified and result.config.sequence.show_track)
        self.seq_check.blockSignals(False)
        self.hits_combo.clear()
        self._search_hits = []

    def _toggle_sequence(self, checked: bool) -> None:
        if self.result is None:
            return
        self.result.config.sequence.show_track = bool(checked)
        self.redraw()

    def _go_to(self) -> None:
        if self.result is None:
            return
        try:
            self.navigator.go_to(self.goto_edit.text())
        except RangeError as e:
            self.cursor_label.setText(str(e))

    def _visible_bases(self) -> tuple[int, int]:
        lo, hi = self.navigator.visible()
        return visible_base_range(lo, hi, self.result.cons_length)

    def _copy_sequence(self) -> None:
        """Copy the visible consensus region as FASTA (explicit action)."""
        from PySide6.QtWidgets import QApplication

        from tescope import __version__

        from ...core.sequence.ops import fasta_text

        if self.result is None or not self._seq_caps().copy_sequence:
            return
        b1, b2 = self._visible_bases()
        text = fasta_text(self.result.sequence.sequence, b1, b2,
                          tescope_version=__version__)
        QApplication.clipboard().setText(text)
        self.cursor_label.setText(
            f"Copied bases {b1:,}-{b2:,} (1-based, {b2 - b1 + 1:,} bp) "
            "as FASTA to the clipboard.")

    def _export_fasta(self) -> None:
        from tescope import __version__

        from ...core.sequence.ops import fasta_text

        if self.result is None or not self._seq_caps().export_fasta:
            return
        seq = self.result.sequence.sequence
        b1, b2 = self._visible_bases()
        full = (b1, b2) == (1, seq.length)
        if not full:
            box = QMessageBox(self)
            box.setWindowTitle("Export FASTA")
            box.setText(f"Export the full consensus ({seq.length:,} bp) or "
                        f"the visible region ({b1:,}-{b2:,})?")
            full_btn = box.addButton("Full consensus", QMessageBox.AcceptRole)
            box.addButton(f"Visible {b1:,}-{b2:,}", QMessageBox.AcceptRole)
            box.addButton(QMessageBox.Cancel)
            box.exec()
            if box.clickedButton() is None or \
                    box.standardButton(box.clickedButton()) == QMessageBox.Cancel:
                return
            if box.clickedButton() is full_btn:
                b1, b2 = 1, seq.length
        path, _ = QFileDialog.getSaveFileName(
            self, "Export consensus FASTA",
            f"{self.result.repeat_name}_consensus.fa", "FASTA (*.fa *.fasta)")
        if not path:
            return
        Path(path).write_text(
            fasta_text(seq, b1, b2, tescope_version=__version__),
            encoding="utf-8")
        self.cursor_label.setText(
            f"Wrote bases {b1:,}-{b2:,} (1-based) to {path}")

    def _search_sequence(self) -> None:
        from ...core.sequence.model import SequenceError
        from ...core.sequence.ops import search_literal

        if self.result is None or not self._seq_caps().search:
            return
        query = self.search_edit.text().strip()
        if not query:
            return
        try:
            hits = search_literal(self.result.sequence.sequence, query)
        except SequenceError as e:
            self.cursor_label.setText(str(e))
            return
        self._search_hits = hits
        self.hits_combo.clear()
        if not hits:
            self.hits_combo.addItem("no hits")
            self.cursor_label.setText(
                f"No literal match for '{query.upper()}' "
                "(forward strand; coordinates are 1-based).")
            return
        self.hits_combo.addItems([h.label for h in hits])
        self.cursor_label.setText(
            f"{len(hits)} hit(s) — select one to zoom "
            "(1-based inclusive coordinates).")
        self._go_to_hit(0)

    def _go_to_hit(self, index: int) -> None:
        if 0 <= index < len(self._search_hits):
            h = self._search_hits[index]
            self.navigator.zoom_to_bases(h.base_start, h.base_end,
                                         margin_frac=2.0)

    def _open_style_dialog(self) -> None:
        if self.result is None:
            QMessageBox.information(
                self, "Figure style", "Run an analysis first — the style "
                "editor works on the current figure's series.")
            return
        from ...plotting.profiles import series_labels
        from .style_dialog import PlotStyleDialog

        def apply_and_redraw():
            self.redraw()
            self._redraw_heatmap()

        conts = self.result.continuous_tracks()
        caps = {
            "has_continuous": bool(conts),
            "has_matrix": any(
                t.matrix is not None and not t.config.group for t in conts),
        }
        dlg = PlotStyleDialog(
            self.result.config.style, series_labels(self.result),
            on_apply=apply_and_redraw, parent=self, input_caps=caps,
        )
        dlg.exec()

    # -- export ------------------------------------------------------------

    def _export_figure(self) -> None:
        if self._current_fig is None or self.result is None:
            return
        path, chosen = QFileDialog.getSaveFileName(
            self, "Export figure", "tescope_profile.png",
            "PNG (*.png);;SVG (*.svg);;PDF (*.pdf)",
        )
        if not path:
            return
        stem, fmt = _resolve_export(path, chosen)
        try:
            save_figure(self._current_fig, stem, [fmt], self.result)
        except (OSError, ValueError) as e:
            QMessageBox.critical(self, "Export failed", str(e))
            return
        QMessageBox.information(
            self, "Export", f"Figure written to {stem}.{fmt}")

    def _export_heatmap(self) -> None:
        if self.heat_canvas.canvas is None or self.result is None:
            return
        path, chosen = QFileDialog.getSaveFileName(
            self, "Export heatmap", "tescope_heatmap.png",
            "PNG (*.png);;SVG (*.svg);;PDF (*.pdf)",
        )
        if not path:
            return
        stem, fmt = _resolve_export(path, chosen)
        try:
            save_figure(self.heat_canvas.canvas.figure, stem, [fmt],
                        self.result)
        except (OSError, ValueError) as e:
            QMessageBox.critical(self, "Export failed", str(e))
            return
        QMessageBox.information(
            self, "Export", f"Heatmap written to {stem}.{fmt}")

    def _export_data(self) -> None:
        if self.result is None:
            return
        directory = QFileDialog.getExistingDirectory(self, "Export result tables to…")
        if not directory:
            return
        files = export_tables(self.result, directory)
        QMessageBox.information(
            self, "Export", f"Wrote {len(files)} files to {directory}"
        )
