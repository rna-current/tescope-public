"""End-to-end golden analysis: config -> orchestrator -> results -> export.

Covers validation milestones A-H on the synthetic dataset, plus a numeric
regression snapshot (tests/golden/expected_profile.tsv) that pins the
published calculation across platforms and future refactors.
"""

from __future__ import annotations

import json
from pathlib import Path

import golden_dataset as gd
import numpy as np
import pandas as pd
import pytest

from tescope.config.io import load_config, save_config
from tescope.config.model import (
    AggregationConfig,
    AnalysisConfig,
    RepeatSelection,
    TrackConfig,
)
from tescope.core.analysis import AnalysisError, run_analysis
from tescope.core.filtering.filters import FilterConfig
from tescope.core.qc import qc_table

BIN = 100
EXPECTED_N_ELEMENTS = [4, 4, 3, 3, 3, 3, 2, 3, 3, 3]


def make_config(golden_dir, **kwargs) -> AnalysisConfig:
    defaults = dict(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=[
            TrackConfig(file=golden_dir["bigwig"], name="signal_bw", group="ctrl", replicate=1),
            TrackConfig(file=golden_dir["bedgraph"], name="signal_bg", group="ctrl", replicate=2),
            TrackConfig(file=golden_dir["peaks_bed"], name="peaks"),
            TrackConfig(file=golden_dir["sites"], name="sites", kind="point"),
        ],
        aggregation=AggregationConfig(bin_size=BIN),
    )
    defaults.update(kwargs)
    return AnalysisConfig(**defaults)


@pytest.fixture(scope="module")
def result(golden_dir, golden_store):
    return run_analysis(make_config(golden_dir), store=golden_store)


class TestMultiTrack:
    def test_all_tracks_processed(self, result):
        # Milestone G: continuous + interval + point tracks render together
        assert len(result.tracks) == 4
        assert len(result.continuous_tracks()) == 2
        assert len(result.interval_tracks()) == 1
        assert len(result.point_tracks()) == 1

    def test_bigwig_and_bedgraph_agree(self, result):
        a, b = (t.profile for t in result.continuous_tracks())
        np.testing.assert_allclose(a["mean"], b["mean"], rtol=1e-6)
        np.testing.assert_allclose(a["median"], b["median"], rtol=1e-6)

    def test_group_combining(self, result):
        assert "ctrl" in result.group_profiles
        g = result.group_profiles["ctrl"]
        assert g["n_replicates"].max() == 2

    def test_profile_values(self, result):
        prof = result.continuous_tracks()[0].profile
        centers = np.arange(10) * BIN + (BIN - 1) / 2
        expected = centers.copy()
        expected[1] = (3 * 149.5 + 124.5) / 4  # SHORT's partial bin
        np.testing.assert_allclose(prof["mean"], expected, rtol=1e-6)
        assert prof["n_elements"].tolist() == EXPECTED_N_ELEMENTS

    def test_peak_projection(self, result):
        proj = result.interval_tracks()[0].projected
        got = {r.peak_name: (r.cons_start, r.cons_end) for r in proj.itertuples()}
        assert got == gd.EXPECTED_PEAK_PROJECTIONS

    def test_point_projection(self, result):
        pts = result.point_tracks()[0].points
        assert sorted(pts["cons_position"]) == pytest.approx(sorted(gd.EXPECTED_SITE_CONS))

    def test_provenance_recorded(self, result):
        assert result.provenance["assembly"] == "testasm"
        assert result.config.annotation_provenance["md5"]

    def test_qc_table(self, result):
        qc = qc_table(result)
        assert {"section", "metric", "value"} == set(qc.columns)
        assert (qc["section"] == "filtering").any()
        peaks_rows = qc[qc["section"] == "track: peaks"]
        assert (peaks_rows["metric"] == "peaks in file").any()


class TestFiltersEndToEnd:
    def test_length_filter(self, golden_dir, golden_store):
        # Milestone D: the 50 bp fragment is excluded
        cfg = make_config(golden_dir, filters=FilterConfig(min_genome_length=60))
        res = run_analysis(cfg, store=golden_store)
        assert res.n_selected == len(gd.COPIES) - 1
        assert 60_000 not in set(res.copies["genome_start"])

    def test_coverage_filter(self, golden_dir, golden_store):
        # Milestone E: low-completeness copies are excluded
        cfg = make_config(golden_dir, filters=FilterConfig(min_consensus_coverage=0.5))
        res = run_analysis(cfg, store=golden_store)
        assert set(res.copies["genome_start"]) == {10_000, 20_000}
        prof = res.continuous_tracks()[0].profile
        assert prof["n_elements"].tolist() == [2] * 10

    def test_no_copies_left_is_clear_error(self, golden_dir, golden_store):
        cfg = make_config(golden_dir, filters=FilterConfig(min_genome_length=10_000))
        with pytest.raises(AnalysisError, match="No TE annotations satisfy"):
            run_analysis(cfg, store=golden_store)

    def test_few_copies_warning(self, golden_dir, golden_store):
        cfg = make_config(golden_dir, filters=FilterConfig(min_consensus_coverage=0.5))
        res = run_analysis(cfg, store=golden_store)
        assert any("copies remain" in w for w in res.warnings)


class TestChromAliasing:
    def test_alias_track_reconciles(self, golden_dir, golden_store):
        # track uses 'T' while the reference uses 'chrT'
        cfg = make_config(
            golden_dir,
            tracks=[TrackConfig(file=golden_dir["bedgraph_alias"], name="alias")],
        )
        res = run_analysis(cfg, store=golden_store)
        prof = res.continuous_tracks()[0].profile
        centers = np.arange(10) * BIN + (BIN - 1) / 2
        assert prof["mean"][0] == pytest.approx((centers[0] * 4) / 4)
        assert any("alias" in w for w in res.warnings)

    def test_alias_disabled_fails_clearly(self, golden_dir, golden_store):
        cfg = make_config(
            golden_dir,
            tracks=[TrackConfig(file=golden_dir["bedgraph_alias"], name="alias")],
            chrom_alias=False,
        )
        with pytest.raises(AnalysisError, match="shares no chromosomes"):
            run_analysis(cfg, store=golden_store)


class TestConfigRoundtrip:
    def test_yaml_roundtrip_reproduces(self, golden_dir, golden_store, tmp_path):
        # Milestone H groundwork: saved config -> identical results
        cfg = make_config(golden_dir)
        p = tmp_path / "analysis.yaml"
        save_config(cfg, p)
        cfg2 = load_config(p)
        res1 = run_analysis(cfg, store=golden_store)
        res2 = run_analysis(cfg2, store=golden_store)
        np.testing.assert_allclose(
            res1.continuous_tracks()[0].profile["mean"],
            res2.continuous_tracks()[0].profile["mean"],
        )

    def test_version_stamped(self, golden_dir, tmp_path):
        from tescope import __version__

        cfg = make_config(golden_dir)
        p = tmp_path / "a.yaml"
        save_config(cfg, p)
        assert f"tescope_version: {__version__}" in p.read_text()


class TestRegressionSnapshot:
    """Pin the golden numbers. If this fails, a calculation changed —
    either fix the regression or consciously regenerate the snapshot
    (scripts/update_golden_snapshot.py) and document why."""

    SNAPSHOT = Path(__file__).parent / "expected_profile.tsv"

    def test_profile_matches_snapshot(self, result):
        prof = result.continuous_tracks()[0].profile
        if not self.SNAPSHOT.exists():
            pytest.skip("snapshot not yet generated")
        exp = pd.read_csv(self.SNAPSHOT, sep="\t")
        for col in ("mean", "median", "sem"):
            np.testing.assert_allclose(
                prof[col], exp[col], rtol=1e-6, equal_nan=True, err_msg=col
            )
        assert prof["n_elements"].tolist() == exp["n_elements"].tolist()
        assert prof["n_observations"].tolist() == exp["n_observations"].tolist()


class TestExports:
    def test_tables_and_figures(self, result, tmp_path):
        from tescope.export.figures import render_result_figures
        from tescope.export.tables import export_tables

        files = export_tables(result, tmp_path)
        names = {f.name for f in files}
        assert "profile_signal_bw.tsv" in names
        assert "occurrence_peaks.tsv" in names
        assert "projected_intervals_peaks.tsv" in names
        assert "points_sites.tsv" in names
        assert "qc.tsv" in names
        assert "copies.tsv" in names
        assert "analysis.yaml" in names
        assert "group_ctrl.tsv" in names

        # profile table carries the spec'd context columns
        prof = pd.read_csv(tmp_path / "profile_signal_bw.tsv", sep="\t")
        assert {"repeat_name", "signal_name", "group", "aggregation",
                "consensus_bin", "consensus_start", "consensus_end",
                "n_elements"} <= set(prof.columns)

        figs = render_result_figures(result, tmp_path)
        fig_names = {f.name for f in figs}
        assert "profile_stacked.png" in fig_names
        assert "profile_stacked.svg" in fig_names
        assert "profile_stacked.figure.json" in fig_names
        # heatmaps rendered for continuous tracks
        assert any(n.startswith("heatmap_") for n in fig_names)

    def test_csv_format(self, result, tmp_path):
        from tescope.export.tables import export_tables

        files = export_tables(result, tmp_path, formats=["csv"])
        assert any(f.suffix == ".csv" for f in files)


class TestAnnotationLayer:
    """Annotation is a metadata/visualization layer: enabling it must not
    change any scientific numerical output (spec: output invariance)."""

    def _annotated_config(self, golden_dir):
        from tescope.config.model import AnnotationConfig

        return make_config(
            golden_dir,
            annotations=AnnotationConfig(
                enabled=True,
                provider="none",  # offline: user features only
                user_features=[
                    {"name": "SYNTH_REGION_1", "start": 100, "end": 450, "label": "demo"},
                    {"name": "SYNTH_REGION_2", "start": 400, "end": 700},
                ],
            ),
        )

    def test_numerical_outputs_invariant(self, golden_dir, golden_store):
        base = run_analysis(make_config(golden_dir), store=golden_store)
        ann = run_analysis(self._annotated_config(golden_dir), store=golden_store)
        for tb, ta in zip(base.continuous_tracks(), ann.continuous_tracks()):
            pd.testing.assert_frame_equal(tb.profile, ta.profile)
        for tb, ta in zip(base.interval_tracks(), ann.interval_tracks()):
            pd.testing.assert_frame_equal(tb.occurrence, ta.occurrence)
            pd.testing.assert_frame_equal(tb.projected, ta.projected)
        for tb, ta in zip(base.point_tracks(), ann.point_tracks()):
            pd.testing.assert_frame_equal(tb.points, ta.points)
        pd.testing.assert_frame_equal(base.copies, ann.copies)
        np.testing.assert_array_equal(base.element_covered, ann.element_covered)
        assert base.filter_report.to_rows() == ann.filter_report.to_rows()

    def test_features_attached_with_qc_and_exports(self, golden_dir, golden_store, tmp_path):
        from tescope.core.qc import qc_table
        from tescope.export.tables import export_tables

        res = run_analysis(self._annotated_config(golden_dir), store=golden_store)
        assert res.features is not None and len(res.features.user_features) == 2
        assert res.features.user_features[0].origin == "user-defined"
        qc = qc_table(res)
        assert (qc["section"] == "annotation").any()
        assert (qc["section"] == "feature: SYNTH_REGION_1").any()
        # coverage QC definition: golden copies covering [100,450)
        row = qc[(qc["section"] == "feature: SYNTH_REGION_1")
                 & (qc["metric"] == "copies fully covering")]["value"].iloc[0]
        # fully covering [100,450): F [0,1000) and R [0,1000) only
        # (P5 [0,400) stops short; SHORT [100,150) is inside the feature)
        assert row.startswith("2 ")
        over = qc[(qc["section"] == "feature: SYNTH_REGION_1")
                  & (qc["metric"] == "copies overlapping")]["value"].iloc[0]
        assert over.startswith("5 ")  # F, R, P5, MID, SHORT (STAR ends at 100)
        files = export_tables(res, tmp_path)
        names = {f.name for f in files}
        assert "features.tsv" in names and "features_used.json" in names
        snap = json.loads((tmp_path / "features_used.json").read_text())
        assert snap["user_features"][0]["name"] == "SYNTH_REGION_1"

    def test_figure_gains_lane_only_when_enabled(self, golden_dir, golden_store):
        from tescope.plotting.profiles import build_figure

        base = run_analysis(make_config(golden_dir), store=golden_store)
        ann = run_analysis(self._annotated_config(golden_dir), store=golden_store)
        fig_base = build_figure(base)
        fig_ann = build_figure(ann)
        assert len(fig_ann.axes) == len(fig_base.axes) + 1
        # hidden when show_in_plots is off
        ann.config.annotations.show_in_plots = False
        fig_hidden = build_figure(ann)
        assert len(fig_hidden.axes) == len(fig_base.axes)

    def test_headless_cli_run_with_annotations(self, golden_dir, golden_store, tmp_path):
        from tescope.config.io import load_config, save_config

        cfg = self._annotated_config(golden_dir)
        p = tmp_path / "annotated.yaml"
        save_config(cfg, p)
        cfg2 = load_config(p)
        res = run_analysis(cfg2, store=golden_store)
        assert res.features is not None and res.features.n_total == 2
