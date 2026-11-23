"""Left-side control panels: genome, TE selection, filters, track manager."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from ...config.model import TrackConfig
from ...core.filtering.filters import COMPLETENESS_CLASSES, FilterConfig
from ...core.tracks.capabilities import probe_capabilities
from ...core.tracks.loader import UnknownTrackFormatError, detect_format
from ..theme import kind_icon

# GUI mapping choices -> (mapping_class, mapping_method); gated per track
# by ACTUAL file capabilities (NH probed from the file, never assumed)
_MAPPING_CHOICES = (
    ("all reads", ("", "nh"), None),
    ("unique (NH=1)", ("unique", "nh"), "nh_tag"),
    ("multi (NH>1)", ("multi", "nh"), "nh_tag"),
    ("unique (MAPQ)", ("unique", "mapq"), "mapq"),
)

_COMPLETENESS_LABELS = {
    "full_length": "Full / near-full length",
    "five_prime_truncated": "5'-truncated",
    "three_prime_truncated": "3'-truncated",
    "internal_fragment": "Internal fragment",
}


class GenomePanel(QGroupBox):
    """Assembly selection + install actions."""

    assembly_changed = Signal(str)
    fetch_requested = Signal(str)  # known assembly to download
    import_requested = Signal()  # custom RepeatMasker import

    def __init__(self, parent=None):
        super().__init__("Reference genome", parent)
        self.setObjectName("panelCard")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(10, 14, 10, 10)
        lay.setSpacing(7)
        row = QHBoxLayout()
        self.combo = QComboBox()
        self.combo.currentTextChanged.connect(self._on_changed)
        row.addWidget(self.combo, 1)
        self.status = QLabel("")
        row.addWidget(self.status)
        lay.addLayout(row)
        btns = QHBoxLayout()
        self.fetch_btn = QPushButton("Download annotation")
        self.fetch_btn.clicked.connect(lambda: self.fetch_requested.emit(self.current()))
        btns.addWidget(self.fetch_btn)
        imp = QPushButton("Import custom…")
        imp.clicked.connect(self.import_requested.emit)
        btns.addWidget(imp)
        lay.addLayout(btns)
        self._installed: set[str] = set()

    def set_assemblies(self, known: list[str], installed: list[str]) -> None:
        self._installed = set(installed)
        current = self.combo.currentText()
        self.combo.blockSignals(True)
        self.combo.clear()
        for name in list(dict.fromkeys(list(known) + list(installed))):
            label = name if name in self._installed else f"{name} (not downloaded)"
            self.combo.addItem(label, userData=name)
        self.combo.blockSignals(False)
        if current:
            idx = self.combo.findData(current)
            if idx >= 0:
                self.combo.setCurrentIndex(idx)
        self._on_changed()

    def current(self) -> str:
        return self.combo.currentData() or ""

    def set_current(self, assembly: str) -> None:
        idx = self.combo.findData(assembly)
        if idx >= 0:
            self.combo.setCurrentIndex(idx)

    def _on_changed(self, *_):
        name = self.current()
        ok = name in self._installed
        self.status.setText("✓ installed" if ok else "not installed")
        self.fetch_btn.setEnabled(not ok and bool(name))
        if name:
            self.assembly_changed.emit(name)


class TePanel(QGroupBox):
    """Class -> family -> name tree with search."""

    repeat_selected = Signal(str)

    def __init__(self, parent=None):
        super().__init__("Transposable element", parent)
        self.setObjectName("panelCard")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(10, 14, 10, 10)
        lay.setSpacing(7)
        self.search = QLineEdit()
        self.search.setPlaceholderText("Search repeat name (e.g. L1HS)…")
        self.search.textChanged.connect(self._refilter)
        lay.addWidget(self.search)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Repeat", "Copies"])
        self.tree.header().setSectionResizeMode(0, QHeaderView.Stretch)
        self.tree.setIndentation(16)
        self.tree.setUniformRowHeights(True)
        self.tree.setMinimumHeight(150)
        self.tree.itemSelectionChanged.connect(self._on_select)
        lay.addWidget(self.tree, 1)
        self.selected_label = QLabel("No repeat selected")
        self.selected_label.setObjectName("hintLabel")
        lay.addWidget(self.selected_label)
        self._hierarchy: pd.DataFrame | None = None

    def set_hierarchy(self, hierarchy: pd.DataFrame | None) -> None:
        self._hierarchy = hierarchy
        self._refilter()

    def current_repeat(self) -> str:
        items = self.tree.selectedItems()
        if items and items[0].childCount() == 0:
            return items[0].data(0, Qt.UserRole) or ""
        return ""

    def select_repeat(self, name: str) -> None:
        matches = self.tree.findItems(name, Qt.MatchExactly | Qt.MatchRecursive, 0)
        for it in matches:
            if it.childCount() == 0:
                self.tree.setCurrentItem(it)
                return

    def _refilter(self, *_):
        self.tree.clear()
        if self._hierarchy is None or len(self._hierarchy) == 0:
            return
        h = self._hierarchy
        q = self.search.text().strip().lower()
        if q:
            h = h[h["repeat_name"].astype(str).str.lower().str.contains(q, regex=False)]
        by_class: dict[str, QTreeWidgetItem] = {}
        by_family: dict[tuple[str, str], QTreeWidgetItem] = {}
        for r in h.itertuples():
            cls, fam, name = str(r.repeat_class), str(r.repeat_family), str(r.repeat_name)
            if cls not in by_class:
                it = QTreeWidgetItem([cls, ""])
                by_class[cls] = it
                self.tree.addTopLevelItem(it)
            fk = (cls, fam)
            if fk not in by_family:
                it = QTreeWidgetItem([fam, ""])
                by_class[cls].addChild(it)
                by_family[fk] = it
            leaf = QTreeWidgetItem([name, f"{int(r.n_copies):,}"])
            leaf.setData(0, Qt.UserRole, name)
            leaf.setTextAlignment(1, Qt.AlignRight | Qt.AlignVCenter)
            leaf.setForeground(1, self.palette().placeholderText())
            by_family[fk].addChild(leaf)
        if q:
            self.tree.expandAll()

    def _on_select(self):
        name = self.current_repeat()
        if name:
            self.selected_label.setText(f"Selected: {name}")
            self.repeat_selected.emit(name)

    def changeEvent(self, event):  # visual only: refresh item colors on theme switch
        super().changeEvent(event)
        if event.type() == event.Type.PaletteChange:
            current = self.current_repeat()
            self._refilter()
            if current:
                self.select_repeat(current)


class FilterPanel(QGroupBox):
    """Fragment filters with live before/after counts."""

    filters_changed = Signal()

    def __init__(self, parent=None):
        super().__init__("TE copy filters", parent)
        self.setObjectName("panelCard")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(10, 14, 10, 10)
        lay.setSpacing(7)

        row = QHBoxLayout()
        self.len_check = QCheckBox("Minimum genomic length (bp)")
        self.len_spin = QSpinBox()
        self.len_spin.setRange(1, 1_000_000)
        self.len_spin.setValue(1000)
        self.len_spin.setEnabled(False)
        self.len_check.toggled.connect(self.len_spin.setEnabled)
        row.addWidget(self.len_check)
        row.addWidget(self.len_spin)
        lay.addLayout(row)

        row = QHBoxLayout()
        self.cov_check = QCheckBox("Minimum consensus coverage")
        self.cov_spin = QDoubleSpinBox()
        self.cov_spin.setRange(0.0, 1.0)
        self.cov_spin.setSingleStep(0.05)
        self.cov_spin.setValue(0.5)
        self.cov_spin.setEnabled(False)
        self.cov_check.toggled.connect(self.cov_spin.setEnabled)
        row.addWidget(self.cov_check)
        row.addWidget(self.cov_spin)
        lay.addLayout(row)

        lay.addWidget(QLabel("Completeness (keep checked classes; none = keep all):"))
        self.comp_checks: dict[str, QCheckBox] = {}
        for cls in COMPLETENESS_CLASSES:
            cb = QCheckBox(_COMPLETENESS_LABELS[cls])
            cb.toggled.connect(lambda *_: self.filters_changed.emit())
            self.comp_checks[cls] = cb
            lay.addWidget(cb)

        for w in (self.len_check, self.cov_check):
            w.toggled.connect(lambda *_: self.filters_changed.emit())
        self.len_spin.valueChanged.connect(lambda *_: self.filters_changed.emit())
        self.cov_spin.valueChanged.connect(lambda *_: self.filters_changed.emit())

        self.counts = QLabel("")
        self.counts.setTextFormat(Qt.RichText)
        lay.addWidget(self.counts)
        self._last_rows: list[tuple[str, int]] | None = None

    def filter_config(self) -> FilterConfig:
        # fields without widgets (include_secondary, strands, end
        # tolerances) round-trip from the last loaded project config
        # instead of silently resetting to defaults
        base = getattr(self, "_loaded_filter_extras", FilterConfig())
        return FilterConfig(
            min_genome_length=self.len_spin.value() if self.len_check.isChecked() else None,
            min_consensus_coverage=(
                self.cov_spin.value() if self.cov_check.isChecked() else None
            ),
            completeness=[c for c, cb in self.comp_checks.items() if cb.isChecked()],
            end_tolerance_frac=base.end_tolerance_frac,
            end_tolerance_bp=base.end_tolerance_bp,
            include_secondary=base.include_secondary,
            strands=list(base.strands or []),
        )

    def apply_filter_config(self, cfg: FilterConfig) -> None:
        self._loaded_filter_extras = cfg
        self.len_check.setChecked(cfg.min_genome_length is not None)
        if cfg.min_genome_length is not None:
            self.len_spin.setValue(cfg.min_genome_length)
        self.cov_check.setChecked(cfg.min_consensus_coverage is not None)
        if cfg.min_consensus_coverage is not None:
            self.cov_spin.setValue(cfg.min_consensus_coverage)
        for c, cb in self.comp_checks.items():
            cb.setChecked(c in (cfg.completeness or []))

    def changeEvent(self, event):  # visual only: recolor counts on theme switch
        super().changeEvent(event)
        if event.type() == event.Type.PaletteChange and self._last_rows is not None:
            self.show_counts(self._last_rows)

    def show_counts(self, rows: list[tuple[str, int]] | None) -> None:
        """Render the live before/after counts with a clear visual hierarchy.

        Same numbers as always — headline shows what will actually be
        analyzed, the per-step breakdown stays fully visible below.
        """
        self._last_rows = rows
        if not rows:
            self.counts.setText("")
            return
        pal = self.palette()
        sec = pal.color(pal.ColorRole.PlaceholderText).name()
        acc = pal.color(pal.ColorRole.Highlight).name()
        total = rows[0][1]
        retained = rows[-1][1]
        if len(rows) == 1:
            html = (
                f"<div style='font-size:14pt; font-weight:600;'>{total:,}</div>"
                f"<div style='color:{sec};'>annotated copies (no filters)</div>"
            )
        else:
            html = (
                f"<div style='font-size:14pt; font-weight:600; color:{acc};'>"
                f"{retained:,}</div>"
                f"<div style='color:{sec};'>retained of {total:,} annotated</div>"
                "<table style='margin-top:6px;'>"
            )
            for label, n in rows[1:]:
                html += (
                    f"<tr><td style='color:{sec}; padding-right:14px;'>{label}</td>"
                    f"<td align='right'>{n:,}</td></tr>"
                )
            html += "</table>"
        self.counts.setText(html)


class TrackPanel(QGroupBox):
    """Track manager: add/remove/rename/group/replicate/show + reorder."""

    tracks_changed = Signal()

    COLS = ("Name", "Kind", "Group", "Rep", "Mapping", "Show")

    def __init__(self, parent=None):
        super().__init__("Signal tracks", parent)
        self.setObjectName("panelCard")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(10, 14, 10, 10)
        lay.setSpacing(7)
        self.table = QTableWidget(0, len(self.COLS))
        self.table.setHorizontalHeaderLabels(self.COLS)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.table.verticalHeader().setVisible(False)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.itemChanged.connect(lambda *_: self.tracks_changed.emit())
        lay.addWidget(self.table, 1)
        hint = QLabel("Drag & drop bigWig / bedGraph / BED / narrowPeak / broadPeak files here")
        hint.setObjectName("hintLabel")
        hint.setWordWrap(True)
        lay.addWidget(hint)
        btns = QHBoxLayout()
        for label, slot in (
            ("Add…", self._add_dialog),
            ("Remove", self._remove_selected),
            ("Up", lambda: self._move(-1)),
            ("Down", lambda: self._move(1)),
        ):
            b = QPushButton(label)
            if label == "Remove":
                b.setObjectName("destructiveButton")
            b.clicked.connect(slot)
            btns.addWidget(b)
        lay.addLayout(btns)
        self._paths: list[str] = []  # row -> file path
        self._extras: list[dict] = []  # row -> no-widget TrackConfig fields
        self._kinds: list[str] = []  # row -> kind string

    # -- public API --------------------------------------------------------

    def add_files(self, paths: list[str]) -> list[str]:
        """Add track files; returns list of error strings for bad files."""
        errors = []
        for p in paths:
            try:
                kind, fmt = detect_format(p)
            except UnknownTrackFormatError as e:
                errors.append(str(e))
                continue
            self._append_row(
                TrackConfig(file=p, name=os.path.basename(p), kind=kind.value, format=fmt)
            )
        self.tracks_changed.emit()
        return errors

    def track_configs(self) -> list[TrackConfig]:
        out = []
        for row in range(self.table.rowCount()):
            rep_text = self.table.item(row, 3).text().strip()
            combo = self.table.cellWidget(row, 4)
            mclass, mmethod = ("", "nh")
            if combo is not None and combo.currentData() is not None:
                mclass, mmethod = combo.currentData()
            extras = (self._extras[row] if row < len(self._extras) else {})
            out.append(
                TrackConfig(
                    file=self._paths[row],
                    name=self.table.item(row, 0).text().strip(),
                    kind=self._kinds[row],
                    group=self.table.item(row, 2).text().strip(),
                    replicate=int(rep_text) if rep_text.isdigit() else None,
                    mapping_class=mclass,
                    mapping_method=mmethod,
                    visible=self.table.item(row, 5).checkState() == Qt.Checked,
                    **extras,
                )
            )
        return out

    def set_track_configs(self, tracks: list[TrackConfig]) -> None:
        self.table.setRowCount(0)
        self._paths, self._kinds = [], []
        self._extras = []  # per-row TrackConfig fields without widgets
        for t in tracks:
            self._append_row(t)
        self.tracks_changed.emit()

    # -- internals ---------------------------------------------------------

    def _append_row(self, t: TrackConfig) -> None:
        self.table.blockSignals(True)
        row = self.table.rowCount()
        self.table.insertRow(row)
        self._paths.append(t.file)
        self._extras.append({
            'sample': t.sample, 'strand': t.strand,
            'color': t.color, 'min_mapq': t.min_mapq,
            'exclude_duplicates': t.exclude_duplicates,
            'unique_mapq': t.unique_mapq, 'format': t.format,
        })
        kind = t.kind or "continuous"
        self._kinds.append(kind)
        name_item = QTableWidgetItem(t.display_name())
        name_item.setToolTip(t.file)
        self.table.setItem(row, 0, name_item)
        kind_item = QTableWidgetItem(kind)
        kind_item.setIcon(kind_icon(kind))
        kind_item.setFlags(kind_item.flags() & ~Qt.ItemIsEditable)
        self.table.setItem(row, 1, kind_item)
        self.table.setItem(row, 2, QTableWidgetItem(t.group))
        self.table.setItem(
            row, 3, QTableWidgetItem("" if t.replicate is None else str(t.replicate))
        )
        self.table.setCellWidget(row, 4, self._mapping_combo(t))
        show = QTableWidgetItem()
        show.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
        show.setCheckState(Qt.Checked if t.visible else Qt.Unchecked)
        self.table.setItem(row, 5, show)
        self.table.blockSignals(False)

    def _mapping_combo(self, t: TrackConfig) -> QComboBox:
        """Capability-gated mapping selector (probes the actual file)."""
        combo = QComboBox()
        try:
            caps = probe_capabilities(t.file, kind=t.kind or None,
                                      fmt=t.format or None)
        except Exception:  # noqa: BLE001 - unknown file: disable everything
            caps = None
        model = combo.model()
        for label, value, need in _MAPPING_CHOICES:
            combo.addItem(label, value)
            item = model.item(combo.count() - 1)
            if need is not None:
                ok = caps is not None and getattr(caps, need, False)
                if not ok:
                    item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
                    item.setToolTip(caps.why_not(need) if caps is not None
                                    else "Track capabilities unknown.")
        if caps is None or not caps.mapping_classification:
            combo.setEnabled(False)
            combo.setToolTip(
                caps.why_not("mapping_classification") if caps is not None
                else "Track capabilities unknown.")
        else:
            combo.setToolTip(
                "Which classified reads contribute coverage (analysis "
                "option, recorded in config/QC). Classification is "
                "aligner-dependent; see docs.")
        # restore saved selection
        target = (t.mapping_class or "", t.mapping_method or "nh")
        for i in range(combo.count()):
            if combo.itemData(i) == target or (
                    target[0] == "" and combo.itemData(i)[0] == ""):
                combo.setCurrentIndex(i)
                if target[0]:
                    break
                if combo.itemData(i)[0] == "":
                    break
        combo.currentIndexChanged.connect(lambda *_: self.tracks_changed.emit())
        return combo

    def _add_dialog(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Add signal tracks",
            "",
            "Signal files (*.bw *.bigWig *.bedGraph *.bg *.bdg *.bed *.narrowPeak "
            "*.broadPeak *.bam *.cram *.gz);;All files (*)",
        )
        if paths:
            for err in self.add_files(paths):
                QMessageBox.warning(self, "Unsupported file", err)

    def _remove_selected(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for r in rows:
            self.table.removeRow(r)
            del self._paths[r]
            del self._kinds[r]
            del self._extras[r]
        self.tracks_changed.emit()  # capability set re-evaluated by consumers

    def _move(self, delta: int):
        row = self.table.currentRow()
        new = row + delta
        if row < 0 or not (0 <= new < self.table.rowCount()):
            return
        cfgs = self.track_configs()
        cfgs[row], cfgs[new] = cfgs[new], cfgs[row]
        self.set_track_configs(cfgs)
        self.table.setCurrentCell(new, 0)


class TeSelectionPanel(QWidget):
    """Composite of genome + TE + filter panels (the left column top)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.genome = GenomePanel()
        self.te = TePanel()
        self.filters = FilterPanel()
        self.annotation = AnnotationPanel()
        self.sequence = SequencePanel()
        lay.addWidget(self.genome)
        lay.addWidget(self.te, 1)
        lay.addWidget(self.filters)
        lay.addWidget(self.annotation)
        lay.addWidget(self.sequence)


class SequencePanel(QGroupBox):
    """Consensus-sequence source controls.

    The sequence dimension is independent of the signal tracks: Automatic
    (Dfam) uses the official API cache-first; user-provided FASTA/library
    files are read locally in place (never copied or redistributed —
    users remain responsible for their libraries' license terms).
    Compatibility with the analysis' consensus model is verified before
    any base is displayed.
    """

    check_requested = Signal(bool)  # refresh flag
    settings_changed = Signal()

    SOURCES = (
        ("Automatic (Dfam)", "dfam"),
        ("User-provided library…", "library"),
        ("User-provided FASTA…", "fasta"),
        ("None", "none"),
    )

    def __init__(self, parent=None):
        super().__init__("Consensus sequence", parent)
        self.setObjectName("panelCard")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(10, 14, 10, 10)
        lay.setSpacing(7)

        self.enable_check = QCheckBox("Attach consensus sequence")
        self.enable_check.toggled.connect(lambda *_: self.settings_changed.emit())
        lay.addWidget(self.enable_check)

        row = QHBoxLayout()
        row.addWidget(QLabel("Source:"))
        self.source_combo = QComboBox()
        for label, _key in self.SOURCES:
            self.source_combo.addItem(label)
        self.source_combo.currentIndexChanged.connect(self._on_source_changed)
        row.addWidget(self.source_combo, 1)
        self.browse_btn = QPushButton("Browse…")
        self.browse_btn.setEnabled(False)
        self.browse_btn.clicked.connect(self._browse)
        row.addWidget(self.browse_btn)
        lay.addLayout(row)

        self.path_label = QLabel("")
        self.path_label.setObjectName("hintLabel")
        self.path_label.setWordWrap(True)
        self.path_label.setVisible(False)
        lay.addWidget(self.path_label)

        srow = QHBoxLayout()
        self.status_label = QLabel("No repeat selected")
        self.status_label.setObjectName("hintLabel")
        self.status_label.setWordWrap(True)
        srow.addWidget(self.status_label, 1)
        self.check_btn = QPushButton("Check source…")
        self.check_btn.setEnabled(False)
        self.check_btn.clicked.connect(lambda: self.check_requested.emit(False))
        srow.addWidget(self.check_btn)
        lay.addLayout(srow)

        self._path = ""

    # -- state --------------------------------------------------------------

    def source_key(self) -> str:
        return self.SOURCES[self.source_combo.currentIndex()][1]

    def set_source_key(self, key: str) -> None:
        for i, (_label, k) in enumerate(self.SOURCES):
            if k == key:
                self.source_combo.setCurrentIndex(i)
                return

    def path(self) -> str:
        return self._path

    def set_path(self, path: str) -> None:
        self._path = path or ""
        self.path_label.setText(self._path)
        self.path_label.setVisible(bool(self._path))

    def set_status(self, text: str, can_check: bool) -> None:
        self.status_label.setText(text)
        self.check_btn.setEnabled(can_check)

    def sequence_config(self):
        from ...config.model import SequenceConfig

        return SequenceConfig(
            enabled=self.enable_check.isChecked(),
            source=self.source_key(),
            path=self._path,
        )

    def apply_sequence_config(self, cfg) -> None:
        self.enable_check.setChecked(cfg.enabled)
        self.set_source_key(cfg.source or "dfam")
        self.set_path(cfg.path)

    # -- internals ----------------------------------------------------------

    def _on_source_changed(self, *_):
        key = self.source_key()
        self.browse_btn.setEnabled(key in ("library", "fasta"))
        self.settings_changed.emit()

    def _browse(self):
        from PySide6.QtWidgets import QFileDialog

        path, _ = QFileDialog.getOpenFileName(
            self, "Consensus FASTA / library",
            "", "FASTA (*.fa *.fasta *.fa.gz *.fasta.gz);;All files (*)",
        )
        if path:
            self.set_path(path)
            self.settings_changed.emit()


class AnnotationPanel(QGroupBox):
    """Consensus-feature annotation controls.

    Database features are looked up per repeat (cache first) and shown
    with source + compatibility state; user-defined features are entered
    in the table or imported from TSV/YAML and are always labeled
    user-defined — never presented as database annotations.
    """

    fetch_requested = Signal(bool)  # refresh flag
    settings_changed = Signal()

    COLS = ("Name", "Start", "End", "Label")

    def __init__(self, parent=None):
        super().__init__("Consensus features", parent)
        self.setObjectName("panelCard")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(10, 14, 10, 10)
        lay.setSpacing(7)

        self.enable_check = QCheckBox("Annotate consensus features")
        self.enable_check.toggled.connect(lambda *_: self.settings_changed.emit())
        lay.addWidget(self.enable_check)

        row = QHBoxLayout()
        self.status_label = QLabel("No repeat selected")
        self.status_label.setObjectName("hintLabel")
        self.status_label.setWordWrap(True)
        row.addWidget(self.status_label, 1)
        self.fetch_btn = QPushButton("Check Dfam…")
        self.fetch_btn.setEnabled(False)
        self.fetch_btn.clicked.connect(lambda: self.fetch_requested.emit(False))
        row.addWidget(self.fetch_btn)
        lay.addLayout(row)

        lay.addWidget(QLabel("User-defined features (0-based, half-open):"))
        self.table = QTableWidget(0, len(self.COLS))
        self.table.setHorizontalHeaderLabels(self.COLS)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        self.table.setMaximumHeight(140)
        self.table.itemChanged.connect(lambda *_: self.settings_changed.emit())
        lay.addWidget(self.table)

        btns = QHBoxLayout()
        for label, slot in (
            ("Add", self._add_row),
            ("Remove", self._remove_selected),
            ("Import…", self._import_dialog),
        ):
            b = QPushButton(label)
            if label == "Remove":
                b.setObjectName("destructiveButton")
            b.clicked.connect(slot)
            btns.addWidget(b)
        self.show_check = QCheckBox("Show in plots")
        self.show_check.setChecked(True)
        self.show_check.toggled.connect(lambda *_: self.settings_changed.emit())
        btns.addWidget(self.show_check)
        lay.addLayout(btns)

        # imported feature FILES (consensus or genome target; validated,
        # previewed, and read in place at Run — docs/feature-import.md)
        frow = QHBoxLayout()
        self.files_label = QLabel("Feature files: none")
        self.files_label.setObjectName("hintLabel")
        self.files_label.setWordWrap(True)
        frow.addWidget(self.files_label, 1)
        add_file_btn = QPushButton("Add feature file…")
        add_file_btn.setToolTip(
            "Import a TE feature/annotation file (tescope-features TSV, "
            "BED, or GFF3) — consensus or genomic coordinates; the file "
            "is validated and previewed before it is attached.")
        add_file_btn.clicked.connect(self._add_feature_file)
        frow.addWidget(add_file_btn)
        clear_files_btn = QPushButton("Clear")
        clear_files_btn.setObjectName("destructiveButton")
        clear_files_btn.clicked.connect(self._clear_feature_files)
        frow.addWidget(clear_files_btn)
        lay.addLayout(frow)
        self._feature_files: list[dict] = []

    # -- feature files -------------------------------------------------------

    def _add_feature_file(self) -> None:
        from PySide6.QtWidgets import QMessageBox

        from ...core.features.files import FeatureFileError, parse_feature_file

        path, _ = QFileDialog.getOpenFileName(
            self, "Add feature file", "",
            "Feature files (*.tsv *.bed *.bed.gz *.gff3 *.gff *.gff.gz "
            "*.gff3.gz);;All files (*)")
        if not path:
            return
        try:
            parsed = parse_feature_file(path)
        except FeatureFileError as e:
            QMessageBox.warning(self, "Feature file invalid", str(e))
            return
        detail = (
            f"Target: {parsed.target}\n"
            f"Reference: {parsed.reference or '(from config/assembly)'}\n"
            f"Coordinates: {parsed.source_convention}\n"
            f"Format: {parsed.format}\n"
            f"Features: {len(parsed.features):,}"
            + ("\nWarnings:\n- " + "\n- ".join(parsed.warnings)
               if parsed.warnings else ""))
        ok = QMessageBox.question(
            self, "Attach feature file?", detail,
            QMessageBox.Yes | QMessageBox.Cancel)
        if ok != QMessageBox.Yes:
            return
        entry: dict = {"path": path}
        if parsed.target:
            entry["target"] = parsed.target
        self._feature_files.append(entry)
        self._refresh_files_label()
        self.settings_changed.emit()

    def _clear_feature_files(self) -> None:
        self._feature_files = []
        self._refresh_files_label()
        self.settings_changed.emit()

    def _refresh_files_label(self) -> None:
        if not self._feature_files:
            self.files_label.setText("Feature files: none")
        else:
            names = [Path(e["path"]).name for e in self._feature_files]
            self.files_label.setText(
                "Feature files: " + ", ".join(names))

    def feature_file_entries(self) -> list[dict]:
        return [dict(e) for e in self._feature_files]

    def set_feature_file_entries(self, entries: list[dict]) -> None:
        self._feature_files = [dict(e) for e in (entries or [])]
        self._refresh_files_label()

    # -- database-side state ------------------------------------------------

    def set_status(self, text: str, can_fetch: bool) -> None:
        self.status_label.setText(text)
        self.fetch_btn.setEnabled(can_fetch)

    # -- user features -------------------------------------------------------

    def _add_row(self, name: str = "", start: str = "", end: str = "",
                 label: str = "") -> None:
        self.table.blockSignals(True)
        r = self.table.rowCount()
        self.table.insertRow(r)
        for c, v in enumerate((name, start, end, label)):
            self.table.setItem(r, c, QTableWidgetItem(str(v)))
        self.table.blockSignals(False)
        self.settings_changed.emit()

    def _remove_selected(self) -> None:
        for r in sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True):
            self.table.removeRow(r)
        self.settings_changed.emit()

    def _import_dialog(self) -> None:
        from ...core.features.user import UserFeatureError, parse_user_features_file

        path, _ = QFileDialog.getOpenFileName(
            self, "Import user features",
            "", "Feature files (*.tsv *.yaml *.yml *.json);;All files (*)",
        )
        if not path:
            return
        try:
            feats = parse_user_features_file(path)
        except UserFeatureError as e:
            from PySide6.QtWidgets import QMessageBox

            QMessageBox.warning(self, "Import failed", str(e))
            return
        for f in feats:
            self._add_row(f.name, f.start, f.end, f.type)

    def user_feature_records(self) -> list[dict]:
        """Rows as config records; malformed rows raise with a row number."""
        from ...core.features.user import UserFeatureError, feature_from_record

        out = []
        for r in range(self.table.rowCount()):
            vals = [
                (self.table.item(r, c).text().strip() if self.table.item(r, c) else "")
                for c in range(len(self.COLS))
            ]
            if not any(vals):
                continue
            rec = {"name": vals[0], "start": vals[1], "end": vals[2]}
            if vals[3]:
                rec["label"] = vals[3]
            try:
                feature_from_record(dict(rec), where=f"row {r + 1}")
            except UserFeatureError as e:
                raise UserFeatureError(f"User feature table: {e}") from e
            rec["start"] = int(rec["start"])
            rec["end"] = int(rec["end"])
            out.append(rec)
        return out

    def set_user_feature_records(self, records: list[dict]) -> None:
        self.table.blockSignals(True)
        self.table.setRowCount(0)
        self.table.blockSignals(False)
        for rec in records or []:
            self._add_row(rec.get("name", ""), rec.get("start", ""),
                          rec.get("end", ""), rec.get("label", ""))
