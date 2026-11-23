"""Offscreen GUI tests for the Loci tab, inspector, and wiring."""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

import golden_dataset as gd  # noqa: E402
from PySide6.QtCore import QItemSelectionModel, Qt  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from tescope.config.model import (  # noqa: E402
    AnalysisConfig,
    RepeatSelection,
    TrackConfig,
)
from tescope.core.analysis import run_analysis  # noqa: E402
from tescope.gui.widgets.inspector import LocusInspector  # noqa: E402
from tescope.gui.widgets.loci import LociTab, LocusTableModel  # noqa: E402
from tescope.gui.widgets.results import ResultsView  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture(scope="module")
def result(golden_dir, golden_store):
    cfg = AnalysisConfig(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file=str(golden_dir["bigwig"]), name="sig"),
                TrackConfig(file=str(golden_dir["narrowpeak"]), name="peaks")],
    )
    cfg.aggregation.bin_size = 100
    return run_analysis(cfg, store=golden_store)


def _select_rows(tab: LociTab, ids: list[str]) -> None:
    sel = tab.table.selectionModel()
    sel.clearSelection()
    for lid in ids:
        row = tab.model.row_of_locus(lid)
        pi = tab.proxy.mapFromSource(tab.model.index(row, 0))
        sel.select(pi, QItemSelectionModel.Select | QItemSelectionModel.Rows)


F_ID = f"{gd.CHROM}:10000-11000(+)"
R_ID = f"{gd.CHROM}:20000-21000(-)"
P3R_ID = f"{gd.CHROM}:40000-40300(-)"


def test_loci_tab_populates_and_counts(qapp, result):
    tab = LociTab()
    tab.set_result(result)
    assert tab.proxy.rowCount() == len(result.copies)
    assert f"{len(result.copies):,}" in tab.status.text()
    assert "no subset applied" in tab.status.text()


def test_model_numeric_sort_and_missing_display(qapp, result):
    model = LocusTableModel(result.locus_table())
    cov_col = next(i for i, (k, _h) in enumerate(model._cols)
                   if k == "consensus_coverage")
    disp = model.data(model.index(0, cov_col), Qt.DisplayRole)
    assert "%" in disp
    raw = model.data(model.index(0, cov_col), Qt.UserRole)
    assert isinstance(raw, float)


def test_filter_hides_rows(qapp, result):
    tab = LociTab()
    tab.set_result(result)
    tab.search_edit.setText("full_length")
    assert 0 < tab.proxy.rowCount() < len(result.copies)
    tab.search_edit.setText("")
    assert tab.proxy.rowCount() == len(result.copies)


def test_selection_changes_nothing_until_applied(qapp, result):
    tab = LociTab()
    tab.set_result(result)
    before = result.continuous_tracks()[0].profile["mean"].copy()
    _select_rows(tab, [F_ID, R_ID])
    assert tab.subset is None
    assert result.config.loci.subset_applied is False
    after = result.continuous_tracks()[0].profile["mean"]
    assert (before.fillna(-1) == after.fillna(-1)).all()


def test_apply_subset_draws_and_records(qapp, result):
    tab = LociTab()
    tab.set_result(result)
    _select_rows(tab, [F_ID, R_ID])
    tab._apply_subset(compare=True)
    assert tab.subset is not None and tab.subset.n_subset == 2
    assert tab.canvas.canvas is not None  # comparison figure drawn
    assert result.config.loci.subset_applied is True
    assert sorted(result.config.loci.selected) == sorted([F_ID, R_ID])
    assert "APPLIED" in tab.status.text()
    tab._clear_subset()
    assert tab.subset is None
    assert result.config.loci.subset_applied is False


def test_collect_and_apply_loci_config_round_trip(qapp, result):
    tab = LociTab()
    tab.set_result(result)
    _select_rows(tab, [F_ID, R_ID])
    tab._apply_subset(compare=False)
    cfg = AnalysisConfig(
        assembly="testasm", repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file="x.bw", name="sig")])
    tab.collect_into(cfg)
    assert cfg.loci.subset_applied is True
    assert sorted(cfg.loci.selected) == sorted([F_ID, R_ID])
    # project-load path (no result yet)
    tab2 = LociTab()
    tab2.apply_loci_config(cfg.loci)
    cfg2 = AnalysisConfig(
        assembly="testasm", repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file="x.bw", name="sig")])
    tab2.collect_into(cfg2)
    assert cfg2.loci.selected == cfg.loci.selected
    assert cfg2.loci.subset_applied is True


def test_stale_persisted_selection_reported(qapp, result, golden_dir,
                                            golden_store):
    cfg = AnalysisConfig(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file=str(golden_dir["bigwig"]), name="sig")],
    )
    cfg.aggregation.bin_size = 100
    cfg.loci.selected = [F_ID, "chrZ:1-2(+)"]
    res = run_analysis(cfg, store=golden_store)
    tab = LociTab()
    recorded = []
    from PySide6.QtWidgets import QMessageBox

    orig = QMessageBox.information
    QMessageBox.information = staticmethod(
        lambda *a, **k: recorded.append(a[2] if len(a) > 2 else ""))
    try:
        tab.set_result(res)
    finally:
        QMessageBox.information = orig
    assert any("chrZ:1-2(+)" in m for m in recorded)
    assert tab.selected_locus_ids() == [F_ID]


def test_inspector_renders_strands_and_projection(qapp, result, golden_dir,
                                                  golden_store, tmp_path):
    cons = tmp_path / "cons.tsv"
    cons.write_text(
        "# tescope-features v1\n# target: consensus\n# reference: TST1\n"
        "feature_id\tstart\tend\tlabel\ttype\n"
        "D1\t100\t400\tdomainX\tdomain\n"
    )
    gen = tmp_path / "genes.bed"
    gen.write_text(f"{gd.CHROM}\t9000\t9500\tneighborGene\t0\t+\n")
    cfg = AnalysisConfig(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file=str(golden_dir["bigwig"]), name="sig")],
    )
    cfg.aggregation.bin_size = 100
    cfg.annotations.feature_files = [
        {"path": str(cons)}, {"path": str(gen), "reference": "testasm"}]
    res = run_analysis(cfg, store=golden_store)
    insp = LocusInspector()
    ids = [F_ID, P3R_ID]
    insp.show_locus(res, ids, F_ID)
    html = insp.body.text()
    assert "10,001–11,000" in html  # 1-based display
    assert "domainX" in html and "projected to" in html
    assert "neighborGene" in html  # within default 5 kb flank
    assert "not gap-aware" in html
    # reverse 5'-truncated copy: domain at cons 100-400 is OUTSIDE its
    # footprint (700-1000) -> no verified correspondence
    insp._step(+1)
    html = insp.body.text()
    assert "no verified correspondence" in html
    assert "(- strand)" in html
    assert insp.pos_label.text().startswith("2 of 2")


def test_highlight_switches_to_heatmap(qapp, result):
    view = ResultsView()
    view.show_result(result)
    view.highlight_locus_in_heatmap(F_ID)
    assert view.tabText(view.currentIndex()) == "Heatmap"
    assert view._heat_highlight is not None


def test_annotation_panel_feature_file_entries_round_trip(qapp):
    from tescope.gui.widgets.panels import AnnotationPanel

    panel = AnnotationPanel()
    assert panel.feature_file_entries() == []
    entries = [{"path": "/x/f.tsv", "target": "consensus"}]
    panel.set_feature_file_entries(entries)
    assert panel.feature_file_entries() == entries
    assert "f.tsv" in panel.files_label.text()
    panel._clear_feature_files()
    assert panel.feature_file_entries() == []


def test_toolbar_save_writes_provenance_sidecar(qapp, result, tmp_path,
                                                monkeypatch):
    view = ResultsView()
    view.show_result(result)
    target = tmp_path / "toolbar_saved.png"
    from PySide6.QtWidgets import QFileDialog

    monkeypatch.setattr(
        QFileDialog, "getSaveFileName",
        staticmethod(lambda *a, **k: (str(target), "PNG (*.png)")))
    view.profile_canvas.toolbar.save_figure()
    assert target.exists()
    sidecar = tmp_path / "toolbar_saved.figure.json"
    assert sidecar.exists(), "toolbar Save must write the provenance sidecar"
    import json

    doc = json.loads(sidecar.read_text())
    assert doc["config"]["repeat"]["name"] == "TST1"
