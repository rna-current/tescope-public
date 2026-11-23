"""Plot-style subsystem tests.

Style is a pure rendering layer: the invariance tests here prove that no
visual option changes any scientific number. Rendering checks assert
artist properties instead of pixels (robust across platforms).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tescope.config.model import (
    AnalysisConfig,
    AnnotationConfig,
    RepeatSelection,
    TrackConfig,
)
from tescope.core.analysis import run_analysis
from tescope.plotting.plotstyle import (
    LINE_STYLES,
    PALETTES,
    PlotStyleConfig,
    PlotStyleError,
)
from tescope.plotting.profiles import build_figure, build_heatmap, series_labels

BIN = 100


def cfg_for(golden_dir, style=None, **kw):
    defaults = dict(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=[
            TrackConfig(file=golden_dir["bigwig"], name="sig_a", group="ctrl", replicate=1),
            TrackConfig(file=golden_dir["bedgraph"], name="sig_b"),
            TrackConfig(file=golden_dir["peaks_bed"], name="peaks"),
        ],
    )
    defaults.update(kw)
    c = AnalysisConfig(**defaults)
    c.aggregation.bin_size = BIN
    if style is not None:
        c.style = style
    return c


HEAVY_STYLE = PlotStyleConfig(
    preset="custom",
    palette="grayscale",
    fig_width_cm=8.5,
    fig_height_cm=12.0,
    title_size=14, label_size=11, tick_size=6, legend_size=10, annotation_size=9,
    title="Custom title", x_label="custom x", y_label="custom y",
    show_legend=True, legend_position="outside right",
    show_grid=True, grid_alpha=0.5,
    y_min=-10.0, y_max=1200.0, x_min=50.0, x_max=900.0,
    line_width=2.5,
    series={"ctrl": {"color": "#123456", "line_style": "dashed", "line_width": 3.0,
                     "alpha": 0.7}},
    band_alpha=0.5,
    annotation_alpha=0.6, annotation_border_width=1.5,
    database_lane_color="#AABBCC", user_lane_color="#CCBBAA",
    heatmap_cmap="magma", heatmap_vmin=0.0, heatmap_vmax=500.0,
    heatmap_colorbar=False,
)


class TestModel:
    def test_defaults_valid_and_roundtrip(self):
        st = PlotStyleConfig()
        st.validate()
        st2 = PlotStyleConfig.from_dict(st.to_dict())
        assert st2 == st

    def test_heavy_roundtrip(self):
        st = PlotStyleConfig.from_dict(HEAVY_STYLE.to_dict())
        assert st.series_color("ctrl", 0) == "#123456"
        assert st.series_line_style("ctrl", 0) == "dashed"
        assert st.series_line_width("ctrl") == 3.0

    @pytest.mark.parametrize("field,value,msg", [
        ("preset", "fancy", "preset"),
        ("palette", "rainbow", "palette"),
        ("legend_position", "under the sofa", "legend_position"),
        ("heatmap_cmap", "jet2", "heatmap_cmap"),
        ("grid_alpha", 1.5, "grid_alpha"),
        ("line_width", -1.0, "line_width"),
        ("database_lane_color", "#GGGGGG", "invalid color"),
    ])
    def test_validation_rejects(self, field, value, msg):
        st = PlotStyleConfig()
        setattr(st, field, value)
        with pytest.raises(PlotStyleError, match=msg):
            st.validate()

    def test_invalid_series_color_rejected(self):
        st = PlotStyleConfig(series={"a": {"color": "not-a-color"}})
        with pytest.raises(PlotStyleError, match="invalid color"):
            st.validate()

    def test_bad_axis_range_rejected(self):
        with pytest.raises(PlotStyleError, match="y_max"):
            PlotStyleConfig(y_min=5.0, y_max=1.0).validate()

    def test_palette_determinism_and_explicit_priority(self):
        st = PlotStyleConfig(palette="colorblind",
                             series={"b": {"color": "#FF0000"}})
        assert st.series_color("a", 0) == PALETTES["colorblind"][0]
        assert st.series_color("b", 1) == "#FF0000"  # explicit never reassigned
        assert st.series_color("c", 2) == PALETTES["colorblind"][2]

    def test_grayscale_cycles_line_styles(self):
        st = PlotStyleConfig(palette="grayscale")
        styles = [st.series_line_style(f"s{i}", i) for i in range(4)]
        assert styles == list(LINE_STYLES)
        # explicit line style wins
        st.series = {"s0": {"line_style": "dotted"}}
        assert st.series_line_style("s0", 0) == "dotted"

    def test_preset_reset(self):
        st = PlotStyleConfig.preset_defaults("publication")
        assert st.preset == "publication"
        assert st.fig_width_in() == pytest.approx(17.8 / 2.54)
        scr = PlotStyleConfig.preset_defaults("screen")
        assert scr.fig_width_in() == pytest.approx(7.2)

    def test_stale_series_entries_are_harmless(self):
        st = PlotStyleConfig(series={"gone-group": {"color": "#112233"}})
        st.validate()
        # unknown label just resolves defaults; stale entry stays stored
        assert st.series_color("present", 0) == PALETTES["default"][0]


class TestConfigIntegration:
    def test_config_without_style_section_loads(self, golden_dir, tmp_path):
        from tescope.config.io import load_config

        p = tmp_path / "old.yaml"
        p.write_text(
            f"assembly: testasm\nrepeat:\n  name: TST1\n"
            f"tracks:\n  - file: {golden_dir['bigwig']}\n"
        )
        cfg = load_config(p)
        assert cfg.style.preset == "screen"  # defaults, renders as before

    def test_style_yaml_roundtrip(self, golden_dir, tmp_path):
        from tescope.config.io import load_config, save_config

        cfg = cfg_for(golden_dir, style=PlotStyleConfig.from_dict(HEAVY_STYLE.to_dict()))
        p = tmp_path / "styled.yaml"
        save_config(cfg, p)
        cfg2 = load_config(p)
        assert cfg2.style.to_dict() == cfg.style.to_dict()

    def test_invalid_style_rejected_at_config_level(self, golden_dir):
        from tescope.config.model import ConfigError

        cfg = cfg_for(golden_dir)
        cfg.style.palette = "rainbow"
        with pytest.raises(ConfigError, match="palette"):
            cfg.validate()


class TestScientificInvariance:
    """Mandatory: no style option may change scientific outputs."""

    def _tables(self, res):
        return (
            [t.profile for t in res.continuous_tracks()],
            [t.occurrence for t in res.interval_tracks()],
            res.copies,
            res.element_covered,
        )

    def test_all_styles_identical_numbers(self, golden_dir, golden_store):
        base = run_analysis(cfg_for(golden_dir), store=golden_store)
        styles = [
            PlotStyleConfig(),  # default
            PlotStyleConfig.preset_defaults("publication"),
            PlotStyleConfig.preset_defaults("screen"),
            PlotStyleConfig.from_dict(HEAVY_STYLE.to_dict()),
        ]
        b_profiles, b_occ, b_copies, b_cov = self._tables(base)
        for st in styles:
            res = run_analysis(cfg_for(golden_dir, style=st), store=golden_store)
            profiles, occ, copies, cov = self._tables(res)
            for a, b in zip(b_profiles, profiles):
                pd.testing.assert_frame_equal(a, b)
            for a, b in zip(b_occ, occ):
                pd.testing.assert_frame_equal(a, b)
            pd.testing.assert_frame_equal(b_copies, copies)
            np.testing.assert_array_equal(b_cov, cov)
            # rendering with the style must not mutate the result either
            build_figure(res, mode="overlay")
            build_figure(res, mode="stacked")
            profiles2, *_ = self._tables(res)
            for a, b in zip(profiles, profiles2):
                pd.testing.assert_frame_equal(a, b)

    def test_display_limits_do_not_touch_tables(self, golden_dir, golden_store):
        st = PlotStyleConfig(y_min=0.0, y_max=10.0, x_min=100.0, x_max=500.0)
        res = run_analysis(cfg_for(golden_dir, style=st), store=golden_store)
        fig = build_figure(res, mode="overlay")
        # view window applied...
        assert fig.axes[-1].get_xlim() == (100.0, 500.0)
        # ...but the full data range is still in the tables
        prof = res.continuous_tracks()[0].profile
        assert prof["consensus_end"].max() == 1000

    def test_heatmap_limits_do_not_rescale_matrix(self, golden_dir, golden_store):
        res = run_analysis(cfg_for(golden_dir), store=golden_store)
        track = res.continuous_tracks()[0]
        before = track.matrix.values.copy()
        res.config.style = PlotStyleConfig(heatmap_cmap="Greys",
                                           heatmap_vmin=0.0, heatmap_vmax=1.0)
        fig = build_heatmap(res, track)
        np.testing.assert_array_equal(track.matrix.values, before)
        im = fig.axes[0].images[0]
        assert im.get_clim() == (0.0, 1.0)
        assert im.get_cmap().name == "Greys"


@pytest.fixture(scope="module")
def render_result(golden_dir, golden_store):
    return run_analysis(cfg_for(golden_dir), store=golden_store)


class TestRendering:
    @pytest.fixture()
    def result(self, render_result):
        return render_result

    def test_series_labels_stable(self, result):
        assert series_labels(result) == ["ctrl", "sig_b"]

    def test_artist_properties_match_style(self, result):
        result.config.style = PlotStyleConfig(series={
            "ctrl": {"color": "#123456", "line_style": "dashed", "line_width": 3.0,
                     "alpha": 0.5},
        })
        fig = build_figure(result, mode="overlay")
        ax = next(a for a in fig.axes if a.lines)
        line = {ln.get_label(): ln for ln in ax.lines}["ctrl"]
        assert line.get_color() == "#123456"
        assert line.get_linestyle() == "--"
        assert line.get_linewidth() == 3.0
        assert line.get_alpha() == 0.5
        other = {ln.get_label(): ln for ln in ax.lines}["sig_b"]
        assert other.get_color() == PALETTES["default"][1]

    def test_legend_positions_render(self, result):
        for pos in ("best", "upper left", "outside right"):
            result.config.style = PlotStyleConfig(legend_position=pos)
            fig = build_figure(result, mode="overlay")
            ax = next(a for a in fig.axes if a.lines)
            assert ax.get_legend() is not None
        result.config.style = PlotStyleConfig(show_legend=False)
        fig = build_figure(result, mode="overlay")
        ax = next(a for a in fig.axes if a.lines)
        assert ax.get_legend() is None

    def test_labels_and_title_override(self, result):
        result.config.style = PlotStyleConfig(
            title="SYNTH title", x_label="SYNTH x", y_label="SYNTH y")
        fig = build_figure(result, mode="overlay")
        assert fig._suptitle.get_text() == "SYNTH title"
        assert fig.axes[-1].get_xlabel() == "SYNTH x"
        assert any(a.get_ylabel() == "SYNTH y" for a in fig.axes)

    def test_figure_size_applied(self, result):
        result.config.style = PlotStyleConfig(fig_width_cm=8.5, fig_height_cm=10.0)
        fig = build_figure(result)
        w, h = fig.get_size_inches()
        assert w == pytest.approx(8.5 / 2.54, abs=1e-6)
        assert h == pytest.approx(10.0 / 2.54, abs=1e-6)

    def test_default_style_matches_classic_geometry(self, result):
        result.config.style = PlotStyleConfig()
        fig = build_figure(result)
        assert fig.get_size_inches()[0] == pytest.approx(7.2)

    def test_annotation_lane_styled_coordinates_untouched(self, golden_dir, golden_store):
        cfg = cfg_for(golden_dir, annotations=AnnotationConfig(
            enabled=True, provider="none",
            user_features=[{"name": "SYNTH_F", "start": 100, "end": 400}]))
        cfg.style = PlotStyleConfig(annotation_alpha=0.5,
                                    user_lane_color="#00FF00",
                                    annotation_border_width=2.0)
        res = run_analysis(cfg, store=golden_store)
        fig = build_figure(res)
        f = res.features.user_features[0]
        assert (f.start, f.end) == (100, 400)  # styling never moves features
        # lane still hatched (origin distinguishable despite custom color)
        lane_ax = fig.axes[0]
        assert any(c.get_hatch() for c in lane_ax.collections)


class TestExport:
    def test_all_formats_export(self, golden_dir, golden_store, tmp_path):
        from tescope.export.figures import save_figure

        res = run_analysis(
            cfg_for(golden_dir, style=PlotStyleConfig(legend_position="outside right")),
            store=golden_store)
        fig = build_figure(res, mode="overlay")
        save_figure(fig, tmp_path / "fig", ["png", "svg", "pdf", "tiff"],
                    res, dpi=150)
        for ext in (".png", ".svg", ".pdf", ".tiff"):
            p = tmp_path / f"fig{ext}"
            assert p.exists() and p.stat().st_size > 0, ext
        assert (tmp_path / "fig.figure.json").exists()

    def test_svg_stays_vector(self, golden_dir, golden_store, tmp_path):
        from tescope.export.figures import save_figure

        res = run_analysis(cfg_for(golden_dir), store=golden_store)
        fig = build_figure(res, mode="overlay")
        save_figure(fig, tmp_path / "v", ["svg"], res)
        svg = (tmp_path / "v.svg").read_text()
        assert "<path" in svg  # vector paths, not an embedded raster
        assert "image/png" not in svg

    def test_invalid_format_rejected(self, golden_dir, golden_store, tmp_path):
        from tescope.export.figures import save_figure

        res = run_analysis(cfg_for(golden_dir), store=golden_store)
        fig = build_figure(res)
        with pytest.raises(ValueError, match="unsupported"):
            save_figure(fig, tmp_path / "x", ["bmp"], res)

    def test_headless_styled_run_end_to_end(self, golden_dir, golden_store, tmp_path):
        """GUI-free reproduction: styled config -> figures + tables."""
        from tescope.export.figures import render_result_figures
        from tescope.export.tables import export_tables

        cfg = cfg_for(golden_dir, style=PlotStyleConfig.from_dict(HEAVY_STYLE.to_dict()))
        cfg.output.figure_formats = ["svg", "pdf"]
        res = run_analysis(cfg, store=golden_store)
        export_tables(res, tmp_path)
        figs = render_result_figures(res, tmp_path)
        assert any(f.suffix == ".svg" for f in figs)
        assert any(f.suffix == ".pdf" for f in figs)
        # style serialized into the exported analysis.yaml
        text = (tmp_path / "analysis.yaml").read_text()
        assert "grayscale" in text and "outside right" in text
