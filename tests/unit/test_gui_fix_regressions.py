"""Regressions for the 2026-10-02 GUI fixes: export path/format
resolution, filter-extras round-trip, and project restore of the
Replicates toggle."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6")

from tescope.config.model import FilterConfig


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app


def test_resolve_export_honors_filter_and_keeps_dotted_stems():
    from tescope.gui.widgets.results import _resolve_export

    stem, fmt = _resolve_export("/t/fig.svg", "PNG (*.png)")
    assert (stem.name, fmt) == ("fig", "svg")  # typed extension wins
    stem, fmt = _resolve_export("/t/myfig", "SVG (*.svg)")
    assert (stem.name, fmt) == ("myfig", "svg")  # filter decides
    stem, fmt = _resolve_export("/t/my.sample.v2", "PDF (*.pdf)")
    assert (stem.name, fmt) == ("my.sample.v2", "pdf")  # dotted stem kept
    stem, fmt = _resolve_export("/t/my.sample.png", "SVG (*.svg)")
    assert (stem.name, fmt) == ("my.sample", "png")


def test_filter_panel_round_trips_non_widget_fields(qapp):
    from tescope.gui.widgets.panels import FilterPanel

    panel = FilterPanel()
    loaded = FilterConfig(min_genome_length=1000, completeness=["full_length"],
                          include_secondary=False, strands=["+"],
                          end_tolerance_frac=0.02, end_tolerance_bp=10)
    panel.apply_filter_config(loaded)
    out = panel.filter_config()
    assert out.include_secondary is False
    assert out.strands == ["+"]
    assert out.end_tolerance_frac == 0.02 and out.end_tolerance_bp == 10
    assert out.min_genome_length == 1000
    assert out.completeness == ["full_length"]


def test_filter_panel_defaults_without_loaded_project(qapp):
    from tescope.gui.widgets.panels import FilterPanel

    out = FilterPanel().filter_config()
    assert out.include_secondary is True and out.strands == []


def test_track_panel_round_trips_non_widget_fields(qapp, golden_dir):
    from tescope.config.model import TrackConfig
    from tescope.gui.widgets.panels import TrackPanel

    panel = TrackPanel()
    loaded = [TrackConfig(file=golden_dir["bigwig"], name="sig",
                          strand="+", min_mapq=10, exclude_duplicates=False,
                          unique_mapq=20, sample="S1", color="#112233")]
    panel.set_track_configs(loaded)
    out = panel.track_configs()[0]
    assert out.strand == "+" and out.min_mapq == 10
    assert out.exclude_duplicates is False and out.unique_mapq == 20
    assert out.sample == "S1" and out.color == "#112233"
