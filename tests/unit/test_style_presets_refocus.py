"""New style presets (2026-09 refocus): publication_minimal,
presentation, high_contrast.

Contract checks: every preset validates, resolves sizes, renders, and
changes NO scientific number; preset_defaults applies the recommended
representation overlay; the screen default remains byte-for-byte the
historical rendering configuration.
"""

from __future__ import annotations

import pandas as pd
import pytest

from tescope.config.model import AnalysisConfig, RepeatSelection, TrackConfig
from tescope.core.analysis import run_analysis
from tescope.plotting.plotstyle import PRESETS, PlotStyleConfig
from tescope.plotting.profiles import build_figure

NEW_PRESETS = ("publication_minimal", "presentation", "high_contrast")


def cfg_for(golden_dir, style=None):
    c = AnalysisConfig(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file=golden_dir["bigwig"], name="sig_a"),
                TrackConfig(file=golden_dir["bedgraph"], name="sig_b")],
    )
    c.aggregation.bin_size = 100
    if style is not None:
        c.style = style
    return c


def test_registry_contains_new_presets():
    for p in NEW_PRESETS:
        assert p in PRESETS


@pytest.mark.parametrize("preset", [p for p in PRESETS])
def test_every_preset_validates_and_resolves(preset):
    st = PlotStyleConfig.preset_defaults(preset)
    st.validate()
    assert st.fig_width_in() > 0
    assert st.resolved_line_width() > 0
    for which in ("title", "label", "tick", "legend", "annotation"):
        assert st.font(which) > 0


def test_reset_overlays_are_display_only_and_applied():
    pres = PlotStyleConfig.preset_defaults("presentation")
    assert pres.representation == "filled" and pres.fill_alpha == 0.18
    minimal = PlotStyleConfig.preset_defaults("publication_minimal")
    assert minimal.representation == "line"
    hc = PlotStyleConfig.preset_defaults("high_contrast")
    assert hc.resolved_line_width() >= 2.0
    # screen stays the exact historical default configuration
    assert (PlotStyleConfig.preset_defaults("screen").to_dict()
            == PlotStyleConfig().to_dict())


def test_presentation_is_visibly_heavier_than_publication():
    pres = PlotStyleConfig.preset_defaults("presentation")
    pub = PlotStyleConfig.preset_defaults("publication")
    assert pres.resolved_line_width() > pub.resolved_line_width()
    assert pres.font("title") > pub.font("title")
    assert pres.fig_width_in() > pub.fig_width_in()


@pytest.mark.parametrize("preset", NEW_PRESETS)
def test_new_presets_change_no_scientific_number(golden_dir, golden_store,
                                                 preset):
    base = run_analysis(cfg_for(golden_dir), store=golden_store)
    res = run_analysis(
        cfg_for(golden_dir, style=PlotStyleConfig.preset_defaults(preset)),
        store=golden_store)
    for a, b in zip(base.continuous_tracks(), res.continuous_tracks()):
        pd.testing.assert_frame_equal(a.profile, b.profile)
    pd.testing.assert_frame_equal(base.copies, res.copies)
    # rendering under the preset must not mutate results either
    build_figure(res, mode="overlay")
    build_figure(res, mode="stacked")
    for a, b in zip(base.continuous_tracks(), res.continuous_tracks()):
        pd.testing.assert_frame_equal(a.profile, b.profile)


@pytest.mark.parametrize("preset", NEW_PRESETS)
def test_new_presets_render_with_expected_weight(golden_dir, golden_store,
                                                 preset):
    st = PlotStyleConfig.preset_defaults(preset)
    res = run_analysis(cfg_for(golden_dir, style=st), store=golden_store)
    fig = build_figure(res, mode="overlay")
    lines = [ln for ax in fig.axes for ln in ax.get_lines()]
    assert lines, "no series rendered"
    widths = {round(ln.get_linewidth(), 2) for ln in lines}
    assert st.resolved_line_width() in widths
