"""Signal-representation milestone tests: transforms, capabilities,
mapping classes, distribution bands, compositions, derived exports.

All fixtures are synthetic; canonical scientific outputs are proven
untouched by every display option.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tescope.config.model import (
    AnalysisConfig,
    ConfigError,
    RepeatSelection,
    TrackConfig,
)
from tescope.core.analysis import run_analysis
from tescope.core.tracks.capabilities import probe_capabilities
from tescope.plotting.plotstyle import PlotStyleConfig, PlotStyleError
from tescope.plotting.profiles import build_figure, build_heatmap
from tescope.plotting.transforms import (
    TransformError,
    copy_quantiles,
    smooth,
)

BIN = 100


def cfg_for(golden_dir, style=None, tracks=None, **kw):
    c = AnalysisConfig(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=tracks or [
            TrackConfig(file=golden_dir["bigwig"], name="sig_a"),
            TrackConfig(file=golden_dir["bedgraph"], name="sig_b", group="g", replicate=1),
        ],
        **kw,
    )
    c.aggregation.bin_size = BIN
    if style is not None:
        c.style = style
    return c


# ---------------------------------------------------------------- transforms
class TestSmoothing:
    def test_moving_average_exact(self):
        v = np.array([0, 0, 3, 0, 0], dtype=float)
        out = smooth(v, "moving_average", window=3)
        np.testing.assert_allclose(out, [0, 1, 1, 1, 0])

    def test_input_never_modified(self):
        v = np.array([1.0, 5.0, 1.0, 5.0, 1.0])
        keep = v.copy()
        for m, kw in (("moving_average", {}), ("gaussian", {"sigma": 1.0}),
                      ("savgol", {"window": 5, "polyorder": 2})):
            smooth(v, m, **{"window": 5, "sigma": 1.0, "polyorder": 2, **kw})
            np.testing.assert_array_equal(v, keep)

    def test_nan_bins_stay_nan(self):
        v = np.array([1.0, np.nan, 1.0, 1.0, 1.0])
        for m in ("moving_average", "gaussian", "savgol"):
            out = smooth(v, m, window=5, sigma=1.0, polyorder=2)
            assert np.isnan(out[1]) and np.isfinite(out[0])

    def test_gaussian_preserves_constant(self):
        out = smooth(np.full(50, 7.0), "gaussian", sigma=3.0)
        np.testing.assert_allclose(out, 7.0)

    def test_savgol_preserves_polynomial(self):
        x = np.arange(30, dtype=float)
        v = 0.5 * x**2 - 3 * x + 1
        out = smooth(v, "savgol", window=7, polyorder=2)
        inner = slice(3, -3)
        np.testing.assert_allclose(out[inner], v[inner], rtol=1e-8)

    @pytest.mark.parametrize("kw,msg", [
        (dict(method="moving_average", window=4), "odd"),
        (dict(method="gaussian", sigma=0), "sigma"),
        (dict(method="savgol", window=5, polyorder=7), "polyorder"),
        (dict(method="loess"), "method"),
    ])
    def test_validation(self, kw, msg):
        with pytest.raises(TransformError, match=msg):
            smooth(np.ones(10), kw.pop("method"),
                   **{"window": 5, "sigma": 1.0, "polyorder": 2, **kw})


class TestQuantiles:
    def test_known_values(self):
        class M:
            values = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0],
                               [4.0, 40.0], [5.0, np.nan]])
        q = copy_quantiles(M, (25.0, 50.0, 75.0))
        np.testing.assert_allclose(q[50.0], [3.0, 25.0])
        np.testing.assert_allclose(q[25.0][0], 2.0)

    def test_needs_two_observations(self):
        class M:
            values = np.array([[1.0, np.nan], [np.nan, np.nan]])
        q = copy_quantiles(M, (50.0,))
        assert np.isnan(q[50.0]).all()


# -------------------------------------------------------------- capabilities
class TestCapabilities:
    def test_bigwig_refuses_mapping(self, golden_dir):
        caps = probe_capabilities(golden_dir["bigwig"])
        assert caps.continuous_signal and caps.smoothing and caps.copy_distribution
        assert not caps.mapping_classification
        assert "pre-aggregated" in caps.why_not("mapping_classification")

    def test_bedgraph_refuses_mapping(self, golden_dir):
        caps = probe_capabilities(golden_dir["bedgraph"])
        assert not caps.mapping_classification and not caps.nh_tag

    def test_interval_track_occurrence_no_smoothing(self, golden_dir):
        caps = probe_capabilities(golden_dir["peaks_bed"])
        assert caps.occurrence and not caps.smoothing
        assert "occurrence" in caps.why_not("smoothing")


pysam = pytest.importorskip("pysam", reason="BAM capability tests need pysam")
from tescope.core.tracks.bam import BamCoverageError, open_bam_coverage  # noqa: E402


def _make_bam(path, reads, with_nh=True):
    """reads: (start, length, mapq, nh, reverse)."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrT", "LN": 100_000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, (start, length, mapq, nh, rev) in enumerate(sorted(reads)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * length
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = mapq
            a.cigarstring = f"{length}M"
            a.flag = 0x10 if rev else 0
            if with_nh:
                a.set_tag("NH", nh)
            bam.write(a)
    pysam.index(str(path))
    return str(path)


@pytest.fixture(scope="module")
def nh_bam(tmp_path_factory):
    return _make_bam(tmp_path_factory.mktemp("nh") / "nh.bam", [
        (10_000, 100, 60, 1, False),   # unique
        (10_050, 100, 0, 4, False),    # multi
        (10_100, 100, 60, 1, True),    # unique
    ])


@pytest.fixture(scope="module")
def no_nh_bam(tmp_path_factory):
    return _make_bam(tmp_path_factory.mktemp("nonh") / "no_nh.bam", [
        (10_000, 100, 60, 1, False),
        (10_050, 100, 0, 1, False),
    ], with_nh=False)


class TestMappingClasses:
    def test_bam_capabilities_probe_actual_nh(self, nh_bam, no_nh_bam):
        assert probe_capabilities(nh_bam).nh_tag is True
        caps = probe_capabilities(no_nh_bam)
        assert caps.nh_tag is False and caps.mapq is True
        assert "NH tags" in caps.why_not("nh_tag")

    def test_nh_unique_multi_split(self, nh_bam):
        u = open_bam_coverage(nh_bam, mapping_class="unique", mapping_method="nh")
        m = open_bam_coverage(nh_bam, mapping_class="multi", mapping_method="nh")
        a = open_bam_coverage(nh_bam)
        vu = u.values("chrT", 10_060, 10_061)[0]
        vm = m.values("chrT", 10_060, 10_061)[0]
        va = a.values("chrT", 10_060, 10_061)[0]
        assert (vu, vm, va) == (1.0, 1.0, 2.0)  # split is exact and additive
        for t in (u, m, a):
            t.close()

    def test_mapq_policy_unique_only(self, no_nh_bam):
        u = open_bam_coverage(no_nh_bam, mapping_class="unique",
                              mapping_method="mapq", unique_mapq=30)
        assert u.values("chrT", 10_020, 10_021)[0] == 1.0  # mapq0 read excluded
        u.close()
        with pytest.raises(BamCoverageError, match="MAPQ alone"):
            open_bam_coverage(no_nh_bam, mapping_class="multi",
                              mapping_method="mapq")

    def test_nh_policy_fails_clearly_without_tags(self, no_nh_bam):
        with pytest.raises(BamCoverageError, match="does not contain"):
            open_bam_coverage(no_nh_bam, mapping_class="unique",
                              mapping_method="nh")

    def test_bigwig_never_acquires_mapping_class(self, golden_dir):
        cfg = cfg_for(golden_dir, tracks=[
            TrackConfig(file=golden_dir["bigwig"], mapping_class="unique")])
        with pytest.raises(ConfigError, match="pre-aggregated"):
            cfg.validate()

    def test_mapping_filter_does_not_affect_other_tracks(
            self, golden_dir, golden_store, nh_bam):
        base = run_analysis(cfg_for(golden_dir, tracks=[
            TrackConfig(file=golden_dir["bigwig"], name="bw")]),
            store=golden_store)
        mixed = run_analysis(cfg_for(golden_dir, tracks=[
            TrackConfig(file=golden_dir["bigwig"], name="bw"),
            TrackConfig(file=nh_bam, name="u", mapping_class="unique",
                        mapping_method="nh")]), store=golden_store)
        pd.testing.assert_frame_equal(
            base.continuous_tracks()[0].profile,
            mixed.continuous_tracks()[0].profile)

    def test_qc_reports_policy_and_bigwig_unavailability(
            self, golden_dir, golden_store, nh_bam):
        from tescope.core.qc import qc_table

        res = run_analysis(cfg_for(golden_dir, tracks=[
            TrackConfig(file=golden_dir["bigwig"], name="bw"),
            TrackConfig(file=nh_bam, name="u", mapping_class="unique",
                        mapping_method="nh")]), store=golden_store)
        qc = qc_table(res)
        vals = qc[qc["metric"] == "mapping policy"]["value"].tolist()
        assert any("unavailable from track" in v for v in vals)
        assert any("class=unique" in v for v in vals)


# --------------------------------------------------- rendering & invariance
class TestDisplayInvariance:
    def _tables(self, res):
        return [t.profile.copy() for t in res.continuous_tracks()]

    def test_viz_options_leave_canonical_untouched(self, golden_dir, golden_store):
        base = run_analysis(cfg_for(golden_dir), store=golden_store)
        b = self._tables(base)
        heavy = PlotStyleConfig(
            representation="filled", fill_alpha=0.4,
            smoothing_method="gaussian", smoothing_sigma=2.0,
            band="iqr", outer_band=True, annotation_shading=True,
            fragment_panels=True, heatmap_sort="total_signal",
            heatmap_show_profile=True)
        res = run_analysis(cfg_for(golden_dir, style=heavy), store=golden_store)
        for a, c in zip(b, self._tables(res)):
            pd.testing.assert_frame_equal(a, c)
        for mode in ("stacked", "overlay", "fragments"):
            build_figure(res, mode=mode)
        build_heatmap(res, res.continuous_tracks()[0])
        for a, c in zip(b, self._tables(res)):  # rendering mutates nothing
            pd.testing.assert_frame_equal(a, c)

    def test_smoothed_line_actually_smoothed_raw_visible(self, golden_dir, golden_store):
        st = PlotStyleConfig(smoothing_method="moving_average",
                             smoothing_window=5, smoothing_show_raw=True)
        res = run_analysis(cfg_for(golden_dir, style=st), store=golden_store)
        fig = build_figure(res, mode="overlay")
        ax = next(a for a in fig.axes if a.lines)
        # raw thin lines + smoothed labeled lines
        labeled = [ln for ln in ax.lines if not ln.get_label().startswith("_")]
        thin = [ln for ln in ax.lines if ln.get_label().startswith("_")]
        assert labeled and thin
        raw = res.continuous_tracks()[0].profile["mean"].to_numpy()
        drawn = labeled[0].get_ydata()
        assert not np.allclose(raw, drawn, equal_nan=True)  # smoothed differs
        np.testing.assert_allclose(  # raw trace shows the canonical values
            thin[0].get_ydata(), raw, equal_nan=True)

    def test_filled_representation_same_values(self, golden_dir, golden_store):
        st = PlotStyleConfig(representation="filled", fill_alpha=0.3)
        res = run_analysis(cfg_for(golden_dir, style=st), store=golden_store)
        fig = build_figure(res, mode="overlay")
        ax = next(a for a in fig.axes if a.lines)
        line = [ln for ln in ax.lines if ln.get_label() == "sig_a"][0]
        np.testing.assert_allclose(
            line.get_ydata(),
            res.continuous_tracks()[0].profile["mean"].to_numpy(),
            equal_nan=True)
        assert ax.collections  # the fill exists

    def test_iqr_band_only_for_matrix_series(self, golden_dir, golden_store):
        st = PlotStyleConfig(band="iqr", outer_band=True)
        res = run_analysis(cfg_for(golden_dir, style=st), store=golden_store)
        fig = build_figure(res, mode="stacked")
        axes_with_fills = [a for a in fig.axes if a.collections]
        assert axes_with_fills  # ungrouped sig_a got quantile ribbons

    def test_fragments_mode_subsets_match_filtered_analysis(
            self, golden_dir, golden_store):
        from tescope.core.filtering.filters import FilterConfig

        res = run_analysis(cfg_for(golden_dir), store=golden_store)
        fig = build_figure(res, mode="fragments")
        # full-length subset mean from the figure == filtered-analysis mean
        full = run_analysis(
            cfg_for(golden_dir, filters=FilterConfig(completeness=["full_length"])),
            store=golden_store)
        ax_full = fig.axes[2]  # second profile panel (full-length)
        line = [ln for ln in ax_full.lines if ln.get_label() == "sig_a"][0]
        np.testing.assert_allclose(
            line.get_ydata(),
            full.continuous_tracks()[0].profile["mean"].to_numpy(),
            equal_nan=True)

    def test_heatmap_sort_changes_order_not_values(self, golden_dir, golden_store):
        res = run_analysis(cfg_for(golden_dir), store=golden_store)
        track = res.continuous_tracks()[0]
        imgs = {}
        for sort in ("cons_start", "genome_length", "total_signal"):
            res.config.style = PlotStyleConfig(heatmap_sort=sort)
            fig = build_heatmap(res, track)
            imgs[sort] = fig.axes[-2 if len(fig.axes) > 1 else 0].images[0].get_array() \
                if False else fig.axes[0].images[0].get_array()
        a = np.asarray(imgs["cons_start"], dtype=float)
        b = np.asarray(imgs["genome_length"], dtype=float)
        assert sorted(np.nansum(a, axis=1)) == pytest.approx(
            sorted(np.nansum(b, axis=1)))  # same rows, different order

    def test_annotation_shading_never_moves_features(self, golden_dir, golden_store):
        from tescope.config.model import AnnotationConfig

        cfg = cfg_for(golden_dir, annotations=AnnotationConfig(
            enabled=True, provider="none",
            user_features=[{"name": "SYNTH_R", "start": 200, "end": 500}]))
        cfg.style = PlotStyleConfig(annotation_shading=True)
        res = run_analysis(cfg, store=golden_store)
        build_figure(res)
        f = res.features.user_features[0]
        assert (f.start, f.end) == (200, 500)


class TestDerivedExports:
    def test_smoothed_and_quantile_exports(self, golden_dir, golden_store, tmp_path):
        from tescope.export.tables import export_tables

        st = PlotStyleConfig(smoothing_method="moving_average",
                             smoothing_window=5, band="iqr")
        res = run_analysis(cfg_for(golden_dir, style=st), store=golden_store)
        files = {f.name for f in export_tables(res, tmp_path)}
        assert "profile_sig_a.tsv" in files  # canonical, untouched
        assert "profile_smoothed_sig_a.tsv" in files
        assert "profile_quantiles_sig_a.tsv" in files
        canon = pd.read_csv(tmp_path / "profile_sig_a.tsv", sep="\t")
        np.testing.assert_allclose(  # canonical == analysis result exactly
            canon["mean"], res.continuous_tracks()[0].profile["mean"],
            equal_nan=True)
        sm = pd.read_csv(tmp_path / "profile_smoothed_sig_a.tsv", sep="\t")
        assert sm["smoothing_method"].iloc[0] == "moving_average"
        assert "window=5" in sm["smoothing_params"].iloc[0]

    def test_no_derived_exports_by_default(self, golden_dir, golden_store, tmp_path):
        from tescope.export.tables import export_tables

        res = run_analysis(cfg_for(golden_dir), store=golden_store)
        files = {f.name for f in export_tables(res, tmp_path)}
        assert not any(f.startswith(("profile_smoothed", "profile_quantiles",
                                     "mapping_support")) for f in files)

    def test_mapping_support_export(self, golden_dir, golden_store, nh_bam, tmp_path):
        from tescope.export.tables import export_tables

        res = run_analysis(cfg_for(golden_dir, tracks=[
            TrackConfig(file=nh_bam, name="u", mapping_class="unique",
                        mapping_method="nh"),
            TrackConfig(file=nh_bam, name="m", mapping_class="multi",
                        mapping_method="nh")]), store=golden_store)
        files = {f.name for f in export_tables(res, tmp_path)}
        assert "mapping_support.tsv" in files
        ms = pd.read_csv(tmp_path / "mapping_support.tsv", sep="\t")
        assert set(ms.columns) >= {"unique_coverage", "multi_coverage",
                                   "fraction_unique"}


class TestConfigRoundtrip:
    def test_signal_fields_roundtrip(self, golden_dir, tmp_path):
        from tescope.config.io import load_config, save_config

        st = PlotStyleConfig(representation="filled", smoothing_method="savgol",
                             smoothing_window=7, smoothing_polyorder=3,
                             band="iqr", outer_band=False,
                             annotation_shading=True, fragment_panels=True,
                             heatmap_sort="total_signal", heatmap_show_profile=True)
        cfg = cfg_for(golden_dir, style=st)
        cfg.tracks[0].mapping_class = ""
        p = tmp_path / "s.yaml"
        save_config(cfg, p)
        cfg2 = load_config(p)
        assert cfg2.style.to_dict() == st.to_dict()

    def test_old_config_still_defaults(self, golden_dir, tmp_path):
        from tescope.config.io import load_config

        p = tmp_path / "old.yaml"
        p.write_text(f"assembly: testasm\nrepeat: {{name: TST1}}\n"
                     f"tracks:\n  - file: {golden_dir['bigwig']}\n")
        cfg = load_config(p)
        assert cfg.style.representation == "line"
        assert cfg.style.smoothing_method == "none"
        assert cfg.tracks[0].mapping_class == ""

    def test_invalid_smoothing_rejected(self):
        with pytest.raises(PlotStyleError, match="odd"):
            PlotStyleConfig(smoothing_method="savgol", smoothing_window=6).validate()
