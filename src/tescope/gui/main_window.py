"""TEscope main window: ties panels, worker threads and results together.

All science goes through :func:`tescope.core.analysis.run_analysis` — the
window only collects an :class:`AnalysisConfig` from the panels and
displays the returned result.
"""

from __future__ import annotations

import copy
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFileDialog,
    QInputDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from tescope import __version__

from ..config.io import load_config, save_config
from ..config.model import AnalysisConfig, AnnotationConfig, RepeatSelection
from ..core.analysis import AnalysisResult
from ..core.annotation.schema import add_derived_columns, drop_inconsistent_hits
from ..core.filtering.filters import apply_filters
from ..core.reference.store import ReferenceError, ReferenceStore
from .theme import DEFAULT_THEME, THEMES, app_icon, apply_theme, repolish
from .widgets.panels import TeSelectionPanel, TrackPanel
from .widgets.results import ResultsView
from .worker import (
    AnalysisWorker,
    AnnotationLookupWorker,
    FetchWorker,
    SequenceLookupWorker,
    start_in_thread,
)


class MainWindow(QMainWindow):
    def __init__(self, store: ReferenceStore | None = None,
                 theme_pref: str = DEFAULT_THEME):
        super().__init__()
        self.setWindowTitle(f"TEscope {__version__}")
        self.resize(1280, 820)
        self.setAcceptDrops(True)
        self.store = store or ReferenceStore()
        self._threads: list = []  # keep worker/thread refs alive
        self._fetch_message = ""
        self._repeat_annotation = None  # derived frame for the selected repeat
        self._theme_pref = theme_pref if theme_pref in THEMES else DEFAULT_THEME

        central = QWidget()
        central_lay = QVBoxLayout(central)
        central_lay.setContentsMargins(0, 0, 0, 0)
        central_lay.setSpacing(0)
        central_lay.addWidget(self._build_header())
        splitter = QSplitter(Qt.Horizontal)
        central_lay.addWidget(splitter, 1)
        self.setCentralWidget(central)

        left = QWidget()
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(10, 10, 6, 10)
        left_lay.setSpacing(10)
        self.selection = TeSelectionPanel()
        left_lay.addWidget(self.selection, 1)
        self.tracks = TrackPanel()
        left_lay.addWidget(self.tracks, 1)
        # Scrollable sidebar: on short windows / high display scaling the
        # dense control column scrolls instead of compressing controls into
        # each other. The primary Run action stays pinned below the scroll
        # area so it is always visible. (Visual robustness only — same
        # widgets, same behavior.)
        from PySide6.QtWidgets import QScrollArea

        left_scroll = QScrollArea()
        left_scroll.setObjectName("leftScroll")
        left_scroll.setWidgetResizable(True)
        left_scroll.setFrameShape(QScrollArea.NoFrame)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        left_scroll.setWidget(left)
        left_scroll.setMinimumWidth(400)

        self.run_btn = QPushButton("Run analysis")
        self.run_btn.setObjectName("primaryButton")
        self.run_btn.setMinimumHeight(36)
        self.run_btn.clicked.connect(self.run_analysis)

        left_col = QWidget()
        left_col_lay = QVBoxLayout(left_col)
        left_col_lay.setContentsMargins(0, 0, 0, 10)
        left_col_lay.setSpacing(8)
        left_col_lay.addWidget(left_scroll, 1)
        run_row = QWidget()
        run_lay = QVBoxLayout(run_row)
        run_lay.setContentsMargins(10, 0, 6, 0)
        run_lay.addWidget(self.run_btn)
        left_col_lay.addWidget(run_row)
        splitter.addWidget(left_col)

        self.results = ResultsView()
        splitter.addWidget(self.results)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([420, 860])

        # menus
        file_menu = self.menuBar().addMenu("&File")
        file_menu.addAction("Open project…", self.open_project)
        file_menu.addAction("Save project…", self.save_project)
        file_menu.addSeparator()
        file_menu.addAction("Quit", self.close)
        self._build_view_menu()
        help_menu = self.menuBar().addMenu("&Help")
        help_menu.addAction("About TEscope", self._about)

        # status bar
        self.status_label = QLabel("Ready")
        self.status_label.setObjectName("statusMessage")
        self.statusBar().addWidget(self.status_label, 1)
        self.progress = QProgressBar()
        self.progress.setMaximumWidth(220)
        self.progress.setVisible(False)
        self.statusBar().addPermanentWidget(self.progress)

        # wiring
        self.selection.genome.assembly_changed.connect(self._on_assembly_changed)
        self.selection.genome.fetch_requested.connect(self._fetch_assembly)
        self.selection.genome.import_requested.connect(self._import_custom)
        self.selection.te.repeat_selected.connect(self._on_repeat_selected)
        self.selection.filters.filters_changed.connect(self._update_filter_counts)
        self.selection.annotation.fetch_requested.connect(self._on_annotation_fetch)
        self.selection.sequence.check_requested.connect(self._on_sequence_check)
        self.selection.sequence.settings_changed.connect(
            self._update_sequence_status)
        self._current_cons_length = 0
        self._plot_style = None  # persists Style-dialog edits across runs
        # project settings that have no GUI control (aggregation, output,
        # chromosome-alias policy, export dpi): kept from the opened
        # project so GUI runs reproduce the file exactly
        self._project_base: AnalysisConfig | None = None

        self._refresh_assemblies()

    # -- chrome (visual only) ----------------------------------------------

    def _build_header(self) -> QWidget:
        from PySide6.QtWidgets import QHBoxLayout

        bar = QWidget()
        bar.setObjectName("headerBar")
        lay = QHBoxLayout(bar)
        lay.setContentsMargins(14, 7, 14, 7)
        lay.setSpacing(10)
        icon_label = QLabel()
        icon_label.setPixmap(app_icon().pixmap(26, 26))
        lay.addWidget(icon_label)
        title = QLabel("TEscope")
        title.setObjectName("headerTitle")
        lay.addWidget(title)
        subtitle = QLabel("Transposable Element Consensus Profiling")
        subtitle.setObjectName("headerSubtitle")
        lay.addWidget(subtitle)
        lay.addStretch(1)
        version = QLabel(f"v{__version__}")
        version.setObjectName("headerSubtitle")
        lay.addWidget(version)
        return bar

    def _build_view_menu(self) -> None:
        from PySide6.QtGui import QActionGroup

        view_menu = self.menuBar().addMenu("&View")
        theme_menu = view_menu.addMenu("Theme")
        group = QActionGroup(self)
        group.setExclusive(True)
        for name in THEMES:
            act = theme_menu.addAction(name.capitalize())
            act.setCheckable(True)
            act.setChecked(name == self._theme_pref)
            act.triggered.connect(lambda _=False, n=name: self._set_theme(n))
            group.addAction(act)

    def _set_theme(self, name: str) -> None:
        from PySide6.QtCore import QSettings
        from PySide6.QtWidgets import QApplication

        self._theme_pref = name
        app = QApplication.instance()
        resolved = apply_theme(app, name)
        app.setWindowIcon(app_icon(dark=resolved == "dark"))
        QSettings("TEscope", "TEscope").setValue("theme", name)

    def _set_status_state(self, state: str) -> None:
        """Tint the status message ('' | running | ok | warn | error). GUI thread only."""
        self.status_label.setProperty("state", state)
        repolish(self.status_label)

    # -- drag & drop -------------------------------------------------------

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        paths = [u.toLocalFile() for u in event.mimeData().urls() if u.isLocalFile()]
        errors = self.tracks.add_files(paths)
        for e in errors:
            QMessageBox.warning(self, "Unsupported file", e)

    # -- assemblies --------------------------------------------------------

    def _refresh_assemblies(self):
        known = list(self.store.known_assemblies())
        installed = self.store.installed_assemblies()
        self.selection.genome.set_assemblies(known, installed)

    def _on_assembly_changed(self, assembly: str):
        self._repeat_annotation = None
        self.selection.filters.show_counts(None)
        if not self.store.is_installed(assembly):
            self.selection.te.set_hierarchy(None)
            self.status_label.setText(
                f"Reference annotation for {assembly} not yet downloaded"
            )
            return
        try:
            hierarchy = self.store.hierarchy(assembly)
        except ReferenceError as e:
            QMessageBox.critical(self, "Annotation error", str(e))
            return
        self.selection.te.set_hierarchy(hierarchy)
        self.status_label.setText(
            f"{assembly}: {int(hierarchy['n_copies'].sum()):,} annotations, "
            f"{len(hierarchy):,} repeat names"
        )

    def _fetch_assembly(self, assembly: str):
        worker = FetchWorker(self.store, assembly)
        self._run_fetch(worker, f"Downloading {assembly} annotation…")

    def _import_custom(self):
        rm_path, _ = QFileDialog.getOpenFileName(
            self, "RepeatMasker annotation",
            "", "RepeatMasker output (*.out *.out.gz *.txt *.txt.gz);;All files (*)",
        )
        if not rm_path:
            return
        name, ok = QInputDialog.getText(
            self, "Assembly name", "Name for this custom assembly:"
        )
        if not ok or not name.strip():
            return
        worker = FetchWorker(self.store, name.strip(), custom_path=rm_path)
        self._run_fetch(worker, f"Importing {Path(rm_path).name}…")

    def _run_fetch(self, worker: FetchWorker, message: str):
        thread = start_in_thread(worker, self)
        self._threads.append((thread, worker))
        self._fetch_message = message
        self.progress.setVisible(True)
        self.progress.setRange(0, 0)
        self.status_label.setText(message)

        # Workers emit from their own thread. Slots must be bound methods of
        # this QObject with explicit QueuedConnection so they execute on the
        # GUI thread — connecting local closures makes PySide call them
        # directly in the worker thread and crashes Qt (observed segfault).
        worker.progress.connect(self._on_fetch_progress, Qt.QueuedConnection)
        worker.finished.connect(self._on_fetch_done, Qt.QueuedConnection)
        worker.failed.connect(self._on_fetch_failed, Qt.QueuedConnection)
        thread.start()

    def _on_fetch_progress(self, stage: str, done: int, total: int):
        self._set_status_state("running")
        if stage == "download" and total:
            self.progress.setRange(0, 100)
            self.progress.setValue(100 * done // total)
        else:
            self.progress.setRange(0, 0)
            self.status_label.setText(f"{self._fetch_message} ({stage})")

    def _on_fetch_done(self, assembly: str):
        self.progress.setVisible(False)
        self._set_status_state("ok")
        self.status_label.setText(f"Installed {assembly}")
        self._refresh_assemblies()
        self.selection.genome.set_current(assembly)

    def _on_fetch_failed(self, msg: str):
        self.progress.setVisible(False)
        self._set_status_state("error")
        self.status_label.setText("Reference installation failed")
        QMessageBox.critical(self, "Reference installation failed", msg)

    # -- TE / filters ------------------------------------------------------

    def _on_repeat_selected(self, repeat_name: str):
        assembly = self.selection.genome.current()
        try:
            sub = self.store.select_repeat(assembly, repeat_name)
        except ReferenceError as e:
            QMessageBox.critical(self, "Selection error", str(e))
            return
        cons_lengths = self.store.consensus_lengths(assembly)
        # same exclusion the analysis applies, so live counts match results
        model_len = int(cons_lengths.get(repeat_name, 0) or 0)
        if model_len > 0:
            sub, _ = drop_inconsistent_hits(sub, model_len)
        self._repeat_annotation = add_derived_columns(sub, cons_lengths)
        self._current_cons_length = model_len
        self._update_filter_counts()
        self._update_annotation_status()
        self._update_sequence_status()

    def _update_filter_counts(self):
        if self._repeat_annotation is None:
            return
        cfg = self.selection.filters.filter_config()
        _, report = apply_filters(self._repeat_annotation, cfg)
        self.selection.filters.show_counts(report.to_rows())

    # -- consensus-feature annotation --------------------------------------

    def _update_annotation_status(self) -> None:
        """Cache-only availability check (fast, no network)."""
        from ..core.features.cache import load_cached
        from ..core.features.model import Compatibility, FeatureProvenance

        panel = self.selection.annotation
        name = self.selection.te.current_repeat()
        if not name or self._current_cons_length <= 0:
            panel.set_status("No repeat selected", can_fetch=False)
            return
        snap = load_cached(name, self._current_cons_length,
                           "dfam", self.store.cache_dir)
        if snap is None:
            panel.set_status(
                "Dfam: not checked yet — features will be looked up on Run "
                "when annotation is enabled.", can_fetch=True)
            return
        pv = FeatureProvenance.from_dict(snap.get("provenance") or {})
        n = len(snap.get("features") or [])
        if pv.compatibility == Compatibility.VERIFIED.value:
            panel.set_status(
                f"Dfam {pv.accession} ({pv.database_version}): verified, "
                f"{n} feature(s), cached.", can_fetch=True)
        elif pv.compatibility == Compatibility.INCOMPATIBLE.value:
            panel.set_status(
                f"Dfam ({pv.database_version}): incompatible consensus — "
                "features withheld.", can_fetch=True)
        else:
            panel.set_status(
                f"Dfam ({pv.database_version}): no matching family.",
                can_fetch=True)

    def _on_annotation_fetch(self, refresh: bool) -> None:
        name = self.selection.te.current_repeat()
        if not name or self._current_cons_length <= 0:
            return
        panel = self.selection.annotation
        panel.set_status("Checking Dfam…", can_fetch=False)
        worker = AnnotationLookupWorker(
            AnnotationConfig(enabled=True, provider="dfam"),
            name, self._current_cons_length, self.store.cache_dir,
            refresh=refresh,
        )
        thread = start_in_thread(worker, self)
        self._threads.append((thread, worker))
        worker.finished.connect(self._on_annotation_done, Qt.QueuedConnection)
        worker.failed.connect(self._on_annotation_failed, Qt.QueuedConnection)
        thread.start()

    def _on_annotation_done(self, fs) -> None:
        self._update_annotation_status()
        self.status_label.setText(f"Annotation: {fs.status}")

    def _on_annotation_failed(self, msg: str) -> None:
        self.selection.annotation.set_status(
            f"Annotation lookup failed: {msg}", can_fetch=True)

    # -- consensus sequence -------------------------------------------------

    def _update_sequence_status(self) -> None:
        """Cache/local-only availability check (fast, no network)."""
        panel = self.selection.sequence
        name = self.selection.te.current_repeat()
        if not name or self._current_cons_length <= 0:
            panel.set_status("No repeat selected", can_check=False)
            return
        seq_cfg = panel.sequence_config()
        if seq_cfg.source == "none":
            panel.set_status("No sequence source.", can_check=False)
            return
        if seq_cfg.source in ("library", "fasta") and not seq_cfg.path:
            panel.set_status("Choose a local FASTA/library file (Browse…).",
                             can_check=False)
            return
        from ..core.sequence.resolve import resolve_sequence

        # cache/local resolution only — the Check button does network work
        res = resolve_sequence(seq_cfg, name, self._current_cons_length,
                               cache_dir=self.store.cache_dir,
                               allow_network=False)
        if res.verified:
            pv = res.sequence.provenance
            panel.set_status(
                f"Verified: {res.detail} — provider {pv.provider}"
                + (f" ({pv.source_version})" if pv.source_version else ""),
                can_check=True)
        elif seq_cfg.source == "dfam" and res.status == "unavailable":
            panel.set_status(
                "Dfam: not checked yet — the sequence will be looked up on "
                "Run, or use Check source now.", can_check=True)
        else:
            panel.set_status(f"{res.status}: {res.detail}", can_check=True)

    def _on_sequence_check(self, refresh: bool) -> None:
        name = self.selection.te.current_repeat()
        if not name or self._current_cons_length <= 0:
            return
        panel = self.selection.sequence
        panel.set_status("Checking sequence source…", can_check=False)
        worker = SequenceLookupWorker(
            panel.sequence_config(), name, self._current_cons_length,
            self.store.cache_dir, refresh=refresh,
        )
        thread = start_in_thread(worker, self)
        self._threads.append((thread, worker))
        worker.finished.connect(self._on_sequence_done, Qt.QueuedConnection)
        worker.failed.connect(self._on_sequence_failed, Qt.QueuedConnection)
        thread.start()

    def _on_sequence_done(self, res) -> None:
        self._update_sequence_status()
        self.status_label.setText(f"Sequence: {res.status} — {res.detail}")

    def _on_sequence_failed(self, msg: str) -> None:
        self.selection.sequence.set_status(
            f"Sequence lookup failed: {msg}", can_check=True)

    # -- project save/load -------------------------------------------------

    def collect_config(self) -> AnalysisConfig:
        cfg = AnalysisConfig(
            assembly=self.selection.genome.current(),
            repeat=RepeatSelection(name=self.selection.te.current_repeat()),
            tracks=self.tracks.track_configs(),
            filters=self.selection.filters.filter_config(),
        )
        if self.results.result is not None:
            cfg.aggregation = self.results.result.config.aggregation
        base = self._project_base
        if base is not None:
            cfg.aggregation = copy.deepcopy(base.aggregation)
            cfg.output = copy.deepcopy(base.output)
            cfg.chrom_alias = base.chrom_alias
            cfg.plot.dpi = base.plot.dpi
        cfg.plot.y_axis = self.results.yaxis_combo.currentText()
        cfg.plot.normalization = self.results.norm_combo.currentText()
        cfg.plot.show_sem = self.results.sem_check.isChecked()
        cfg.plot.show_replicates = self.results.reps_check.isChecked()
        if self._plot_style is not None:
            cfg.style = self._plot_style
        self._plot_style = cfg.style  # dialog edits mutate this in place
        panel = self.selection.annotation
        cfg.annotations = AnnotationConfig(
            enabled=panel.enable_check.isChecked(),
            provider="dfam",
            user_features=panel.user_feature_records(),
            feature_files=panel.feature_file_entries(),
            show_in_plots=panel.show_check.isChecked(),
        )
        cfg.sequence = self.selection.sequence.sequence_config()
        self.results.loci_tab.collect_into(cfg)
        if self.results.result is not None:
            # preserve the browser's Show-sequence toggle across runs
            cfg.sequence.show_track = \
                self.results.result.config.sequence.show_track
        return cfg

    def apply_config(self, cfg: AnalysisConfig) -> None:
        self._project_base = cfg
        self.selection.genome.set_current(cfg.assembly)
        self.selection.te.select_repeat(cfg.repeat.name)
        self.selection.filters.apply_filter_config(cfg.filters)
        self.tracks.set_track_configs(cfg.tracks)
        self.results.yaxis_combo.setCurrentText(cfg.plot.y_axis)
        self.results.norm_combo.setCurrentText(cfg.plot.normalization)
        self.results.sem_check.setChecked(cfg.plot.show_sem)
        self.results.reps_check.setChecked(cfg.plot.show_replicates)
        self._plot_style = cfg.style
        panel = self.selection.annotation
        panel.enable_check.setChecked(cfg.annotations.enabled)
        panel.show_check.setChecked(cfg.annotations.show_in_plots)
        panel.set_user_feature_records(cfg.annotations.user_features)
        panel.set_feature_file_entries(cfg.annotations.feature_files)
        self.selection.sequence.apply_sequence_config(cfg.sequence)
        self.results.loci_tab.apply_loci_config(cfg.loci)

    def save_project(self):
        try:
            cfg = self.collect_config()
            cfg.validate()
        except Exception as e:  # noqa: BLE001 - validation message to user
            QMessageBox.warning(self, "Incomplete project", str(e))
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save TEscope project", "analysis.yaml",
            "TEscope project (*.yaml *.yml *.json)",
        )
        if path:
            save_config(cfg, path)
            self.status_label.setText(f"Project saved: {path}")

    def open_project(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open TEscope project", "", "TEscope project (*.yaml *.yml *.json)"
        )
        if not path:
            return
        try:
            cfg = load_config(path)
        except Exception as e:  # noqa: BLE001
            QMessageBox.critical(self, "Cannot open project", str(e))
            return
        self.apply_config(cfg)
        self.status_label.setText(f"Project loaded: {path}")

    # -- analysis ----------------------------------------------------------

    def run_analysis(self):
        try:
            cfg = self.collect_config()
            cfg.validate()
        except Exception as e:  # noqa: BLE001
            QMessageBox.warning(self, "Cannot run", str(e))
            return
        if not self.store.is_installed(cfg.assembly):
            QMessageBox.warning(
                self, "Cannot run",
                f"Reference annotation for {cfg.assembly} is not yet downloaded.",
            )
            return
        self.run_btn.setEnabled(False)
        self.progress.setVisible(True)
        worker = AnalysisWorker(cfg, self.store)
        thread = start_in_thread(worker, self)
        self._threads.append((thread, worker))
        # bound-method slots + QueuedConnection: see _run_fetch for why
        worker.progress.connect(self._on_analysis_progress, Qt.QueuedConnection)
        worker.finished.connect(self._on_result, Qt.QueuedConnection)
        worker.failed.connect(self._on_failed, Qt.QueuedConnection)
        thread.start()

    def _on_analysis_progress(self, msg: str, step: int, total: int):
        self._set_status_state("running")
        self.progress.setRange(0, total)
        self.progress.setValue(step)
        self.status_label.setText(msg)

    def _on_result(self, result: AnalysisResult):
        self.run_btn.setEnabled(True)
        self.progress.setVisible(False)
        self.results.show_result(result)
        msg = (
            f"{result.repeat_name}: {result.n_selected:,}/{result.n_annotated:,} "
            "copies analyzed"
        )
        if result.warnings:
            # non-modal by design: notes are informational and already appear
            # in the QC tab; a blocking dialog per run is unnecessary friction
            msg += f"  —  {len(result.warnings)} note(s), see QC tab"
            self._set_status_state("warn")
        else:
            self._set_status_state("ok")
        self.status_label.setText(msg)

    def _on_failed(self, message: str):
        self.run_btn.setEnabled(True)
        self.progress.setVisible(False)
        self._set_status_state("error")
        self.status_label.setText("Analysis failed")
        QMessageBox.critical(self, "Analysis failed", message)

    def _about(self):
        from .widgets.about_dialog import AboutDialog

        AboutDialog(self).exec()
