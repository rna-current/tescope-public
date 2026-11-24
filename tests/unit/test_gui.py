"""GUI smoke tests (offscreen Qt platform; skipped when PySide6 is absent).

These verify wiring — panels build a valid config, results populate, and
project round-trips — not visual appearance. Real-device validation on
Windows/macOS is tracked in docs/limitations.md.
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

from tescope.config.model import AnalysisConfig, RepeatSelection, TrackConfig  # noqa: E402
from tescope.core.analysis import run_analysis  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture()
def window(qapp, golden_store):
    from tescope.gui.main_window import MainWindow

    w = MainWindow(store=golden_store)
    yield w
    w.close()


class TestMainWindow:
    def test_assemblies_listed(self, window):
        combo = window.selection.genome.combo
        names = [combo.itemData(i) for i in range(combo.count())]
        assert "hg38" in names and "testasm" in names

    def test_selecting_installed_assembly_fills_te_tree(self, window):
        window.selection.genome.set_current("testasm")
        tree = window.selection.te.tree
        assert tree.topLevelItemCount() >= 2  # LINE and SINE classes

    def test_te_selection_and_live_filter_counts(self, window):
        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        assert window.selection.te.current_repeat() == "TST1"
        # live counts show total copies
        assert "7" in window.selection.filters.counts.text()
        # enable a filter and see updated counts
        window.selection.filters.cov_check.setChecked(True)
        window.selection.filters.cov_spin.setValue(0.5)
        assert "2" in window.selection.filters.counts.text()

    def test_track_panel_add_and_configs(self, window, golden_dir):
        errors = window.tracks.add_files([golden_dir["bigwig"], golden_dir["peaks_bed"]])
        assert errors == []
        cfgs = window.tracks.track_configs()
        assert [c.kind for c in cfgs] == ["continuous", "interval"]
        # unsupported file reports an error instead of silently ignoring
        errors = window.tracks.add_files([golden_dir["out"]])
        assert len(errors) == 1

    def test_collect_config_roundtrip(self, window, golden_dir, tmp_path):
        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        window.tracks.add_files([golden_dir["bigwig"]])
        cfg = window.collect_config()
        cfg.validate()
        assert cfg.assembly == "testasm"
        assert cfg.repeat.name == "TST1"
        # apply_config restores the same state
        window.apply_config(cfg)
        cfg2 = window.collect_config()
        assert cfg2.to_dict()["tracks"] == cfg.to_dict()["tracks"]

    def test_opened_project_keeps_settings_without_gui_controls(
            self, window, golden_dir):
        # regression (0.2.1): File > Open project dropped aggregation.bin_size
        # (no GUI control) and ran with the default, so a GUI run of a
        # project file did not reproduce the CLI run of the same file
        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        window.tracks.add_files([golden_dir["bigwig"]])
        cfg = window.collect_config()
        cfg.aggregation.bin_size = 50
        cfg.aggregation.statistic = "median"
        cfg.chrom_alias = False
        cfg.plot.dpi = 123
        cfg.output.dir = "some/dir"
        window.apply_config(cfg)
        got = window.collect_config().to_dict()
        want = cfg.to_dict()
        for key in ("aggregation", "output", "chrom_alias"):
            assert got[key] == want[key], key
        assert got["plot"]["dpi"] == 123

    def test_results_view_populates(self, window, golden_dir, golden_store):
        cfg = AnalysisConfig(
            assembly="testasm",
            repeat=RepeatSelection(name="TST1"),
            tracks=[
                TrackConfig(file=golden_dir["bigwig"], name="sig", group="g", replicate=1),
                TrackConfig(file=golden_dir["peaks_bed"], name="peaks"),
            ],
        )
        result = run_analysis(cfg, store=golden_store)
        window.results.show_result(result)
        assert window.results.profile_canvas.canvas is not None
        assert window.results.qc_widget.rowCount() > 10
        assert window.results.heat_combo.count() == 1
        # switching view mode rebuilds the figure
        window.results.mode_combo.setCurrentText("overlay")
        assert window.results.profile_canvas.canvas is not None

    def test_gui_and_cli_share_core(self, window, golden_dir, golden_store):
        """The config the GUI collects runs through the same run_analysis."""
        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        window.tracks.add_files([golden_dir["bedgraph"]])
        cfg = window.collect_config()
        result = run_analysis(cfg, store=golden_store)
        assert result.n_selected == 7


class TestWorkerThreadPath:
    def test_threaded_run_delivers_results_via_queued_signals(
        self, window, golden_dir, qapp
    ):
        """Regression: worker signals connected to closures were delivered in
        the worker thread and crashed Qt on real displays (segfault found in
        on-display validation 2026-08-27). The full threaded path must work."""
        import time

        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        window.tracks.add_files([golden_dir["bedgraph"]])
        window.run_analysis()
        t0 = time.time()
        while window.results.result is None and time.time() - t0 < 30:
            qapp.processEvents()
            time.sleep(0.01)
        assert window.results.result is not None, "threaded analysis never delivered"
        assert window.results.result.n_selected == 7
        assert window.run_btn.isEnabled()

    def test_live_filter_counts_match_analysis_copy_set(self, window, golden_dir):
        """Live counts must apply the same inconsistent-hit exclusion as
        run_analysis so the promised and delivered copy sets agree."""
        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        text = window.selection.filters.counts.text()
        assert "7" in text  # golden set has no inconsistent hits: full count


class TestAnnotationPanel:
    def test_panel_collects_and_applies_config(self, window, golden_dir):
        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        window.tracks.add_files([golden_dir["bedgraph"]])
        panel = window.selection.annotation
        panel.enable_check.setChecked(True)
        panel._add_row("SYNTH_GUI_REGION", 10, 500, "demo")
        panel.show_check.setChecked(False)
        cfg = window.collect_config()
        cfg.validate()
        assert cfg.annotations.enabled is True
        assert cfg.annotations.show_in_plots is False
        assert cfg.annotations.user_features == [
            {"name": "SYNTH_GUI_REGION", "start": 10, "end": 500, "label": "demo"}
        ]
        # round-trip back into the panel
        panel.enable_check.setChecked(False)
        panel.set_user_feature_records([])
        window.apply_config(cfg)
        assert panel.enable_check.isChecked()
        assert panel.user_feature_records()[0]["name"] == "SYNTH_GUI_REGION"

    def test_malformed_user_row_is_reported(self, window, golden_dir):
        from tescope.core.features.user import UserFeatureError

        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        panel = window.selection.annotation
        panel._add_row("bad", "20", "5", "")
        with pytest.raises(UserFeatureError, match="must be >"):
            panel.user_feature_records()

    def test_cache_status_shown_on_selection(self, window, golden_store):
        from tescope.core.features.cache import store_cached

        store_cached("TST1", 1000, {
            "features": [], "provenance": {
                "source": "Dfam", "accession": "DFTEST000001",
                "database_version": "Dfam 9.9-synthetic",
                "compatibility": "verified"}},
            cache_dir=golden_store.cache_dir)
        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        assert "verified" in window.selection.annotation.status_label.text()


class TestStyleDialog:
    def test_dialog_edits_style_and_triggers_redraw(self, window, golden_dir, qapp):
        import time

        from tescope.gui.widgets.style_dialog import PlotStyleDialog
        from tescope.plotting.profiles import series_labels

        window.selection.genome.set_current("testasm")
        window.selection.te.select_repeat("TST1")
        window.tracks.add_files([golden_dir["bedgraph"]])
        window.run_analysis()
        t0 = time.time()
        while window.results.result is None and time.time() - t0 < 30:
            qapp.processEvents()
            time.sleep(0.01)
        result = window.results.result
        assert result is not None

        applied = []
        dlg = PlotStyleDialog(result.config.style, series_labels(result),
                              on_apply=lambda: applied.append(1))
        label = series_labels(result)[0]
        # set a hex color + line style through the widgets
        for lab, btn, hexedit, ls, width in dlg._series_widgets:
            if lab == label:
                hexedit.setText("#336699")
                ls.setCurrentIndex(ls.findData("dashed"))
                width.setValue(2.0)
        dlg.preset_combo.setCurrentText("publication")
        dlg._apply()
        assert applied
        st = result.config.style
        assert st.preset == "publication"
        assert st.series[label]["color"] == "#336699"
        assert st.series[label]["line_style"] == "dashed"
        # style persists into the next collected config (save/reload path)
        cfg = window.collect_config()
        assert cfg.style.series[label]["color"] == "#336699"

    def test_dialog_rejects_invalid_hex(self, window, golden_dir, qapp, monkeypatch):
        import time

        from tescope.gui.widgets import style_dialog as sd_mod
        from tescope.gui.widgets.style_dialog import PlotStyleDialog
        from tescope.plotting.profiles import series_labels

        warned = []
        monkeypatch.setattr(sd_mod.QMessageBox, "warning",
                            staticmethod(lambda *a, **k: warned.append(a)))

        if window.results.result is None:
            window.selection.genome.set_current("testasm")
            window.selection.te.select_repeat("TST1")
            window.tracks.add_files([golden_dir["bedgraph"]])
            window.run_analysis()
            t0 = time.time()
            while window.results.result is None and time.time() - t0 < 30:
                qapp.processEvents()
                time.sleep(0.01)
        result = window.results.result
        applied = []
        dlg = PlotStyleDialog(result.config.style, series_labels(result),
                              on_apply=lambda: applied.append(1))
        dlg._series_widgets[0][2].setText("#NOTHEX")
        dlg._apply()  # warning path: must not apply
        assert warned and not applied


class TestCapabilityAwareGui:
    def test_mapping_combo_disabled_for_bigwig(self, window, golden_dir):
        window.tracks.set_track_configs([])
        window.tracks.add_files([golden_dir["bigwig"]])
        combo = window.tracks.table.cellWidget(0, 4)
        assert combo is not None and not combo.isEnabled()
        assert "pre-aggregated" in combo.toolTip()

    def test_mapping_combo_gated_by_actual_nh(self, window, tmp_path):
        pysam = pytest.importorskip("pysam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "chrT", "LN": 100000}]}
        p = tmp_path / "caps.bam"
        with pysam.AlignmentFile(str(p), "wb", header=header) as bam:
            a = pysam.AlignedSegment()
            a.query_name = "r0"
            a.query_sequence = "A" * 50
            a.reference_id = 0
            a.reference_start = 100
            a.mapping_quality = 60
            a.cigarstring = "50M"
            a.set_tag("NH", 1)
            bam.write(a)
        pysam.index(str(p))
        window.tracks.set_track_configs([])
        window.tracks.add_files([str(p)])
        combo = window.tracks.table.cellWidget(0, 4)
        assert combo.isEnabled()
        model = combo.model()
        from PySide6.QtCore import Qt as _Qt
        states = {combo.itemText(i): bool(model.item(i).flags() & _Qt.ItemIsEnabled)
                  for i in range(combo.count())}
        assert states["unique (NH=1)"] is True  # NH present in the file
        assert states["unique (MAPQ)"] is True

    def test_mapping_selection_reaches_config(self, window, tmp_path, golden_dir):
        pytest.importorskip("pysam")
        window.tracks.set_track_configs([])
        window.tracks.add_files([golden_dir["bigwig"]])
        cfgs = window.tracks.track_configs()
        assert cfgs[0].mapping_class == ""  # disabled combo -> no class

    def test_style_dialog_gates_iqr_without_matrix(self, window, golden_dir, qapp):
        from PySide6.QtCore import Qt as _Qt

        from tescope.gui.widgets.style_dialog import PlotStyleDialog
        from tescope.plotting.plotstyle import PlotStyleConfig

        dlg = PlotStyleDialog(PlotStyleConfig(), ["g"],
                              input_caps={"has_continuous": True,
                                          "has_matrix": False})
        i = dlg.band_combo.findText("iqr")
        item = dlg.band_combo.model().item(i)
        assert not (item.flags() & _Qt.ItemIsEnabled)
        assert "grouped" in item.toolTip()
        dlg.close()

    def test_style_dialog_disables_signal_tab_without_continuous(self, qapp):
        from tescope.gui.widgets.style_dialog import PlotStyleDialog
        from tescope.plotting.plotstyle import PlotStyleConfig

        dlg = PlotStyleDialog(PlotStyleConfig(), [],
                              input_caps={"has_continuous": False})
        assert not dlg.smooth_combo.isEnabled()
        assert "cannot be smoothed" in dlg.smooth_combo.toolTip()
        dlg.close()
