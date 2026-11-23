"""Arial-first figure typography (2026-09 policy).

The declared family chain is Arial -> Liberation Sans -> DejaVu Sans
-> installed CJK fallbacks -> sans-serif, identical for every style
preset; SVG exports carry the declaration as font-family metadata
(svg.fonttype 'none'); PDF embeds the resolved face (fonttype 42).
"""

from __future__ import annotations

import pytest

from tescope.plotting.plotstyle import PRESETS, PlotStyleConfig
from tescope.plotting.style import (
    PREFERRED_FAMILIES,
    STYLE,
    font_family_chain,
)


def test_declared_chain_is_arial_first():
    chain = font_family_chain()
    assert chain[:3] == ["Arial", "Liberation Sans", "DejaVu Sans"]
    assert chain[-1] == "sans-serif"
    assert STYLE["font.family"] == chain
    assert PREFERRED_FAMILIES[0] == "Arial"


def test_svg_and_pdf_backend_settings_preserve_fonts():
    assert STYLE["svg.fonttype"] == "none"  # text stays text w/ metadata
    assert STYLE["pdf.fonttype"] == 42  # TrueType embedding


@pytest.mark.parametrize("preset", PRESETS)
def test_no_preset_overrides_the_family(preset):
    rc = PlotStyleConfig.preset_defaults(preset).rc_overrides()
    assert "font.family" not in rc
    merged = {**STYLE, **rc}
    assert merged["font.family"][0] == "Arial"


def test_rendered_svg_declares_arial_first(golden_dir, golden_store,
                                           tmp_path):
    from tescope.config.model import (
        AnalysisConfig,
        RepeatSelection,
        TrackConfig,
    )
    from tescope.core.analysis import run_analysis
    from tescope.export.figures import save_figure
    from tescope.plotting.profiles import build_figure

    cfg = AnalysisConfig(
        assembly="testasm", repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file=golden_dir["bigwig"], name="sig_a")])
    cfg.aggregation.bin_size = 100
    res = run_analysis(cfg, store=golden_store)
    fig = build_figure(res, mode="overlay")
    out = tmp_path / "fig"
    save_figure(fig, out, formats=["svg", "pdf", "png"])
    svg = (tmp_path / "fig.svg").read_text()
    assert "Arial" in svg and "font-family" in svg
    first_decl = svg.split("font-family:")[1][:80]
    assert first_decl.strip().strip("'\"").startswith("Arial")
    assert (tmp_path / "fig.pdf").stat().st_size > 0
    assert (tmp_path / "fig.png").stat().st_size > 0
