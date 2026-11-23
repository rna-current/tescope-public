"""Signal projection onto consensus bins (validation milestones A, B, C, F)."""

from __future__ import annotations

import golden_dataset as gd
import numpy as np
import pytest

from tescope.core.annotation.repeatmasker import parse_repeatmasker_out
from tescope.core.annotation.schema import add_derived_columns
from tescope.core.coordinates.mapping import make_bins
from tescope.core.projection.continuous import project_continuous
from tescope.core.projection.coverage import element_coverage
from tescope.core.projection.intervals import (
    interval_occurrence,
    project_intervals,
    project_points,
)
from tescope.core.tracks.bed import open_intervals, open_points
from tescope.core.tracks.bedgraph import open_bedgraph

BIN = 100


@pytest.fixture(scope="module")
def copies(golden_dir):
    df = parse_repeatmasker_out(golden_dir["out"])
    df = df[df["repeat_name"] == "TST1"]
    return add_derived_columns(df)


def _copy_row(copies, gs):
    return copies[copies["genome_start"] == gs]


def bin_centers():
    return np.arange(10) * BIN + (BIN - 1) / 2


class TestContinuousProjection:
    def test_forward_full_copy(self, golden_dir, copies):
        # Milestone A: forward strand copy projects the ramp exactly
        track = open_bedgraph(golden_dir["bedgraph"])
        m = project_continuous(track, _copy_row(copies, 10_000), gd.CONS_LEN, BIN)
        assert m.values.shape == (1, 10)
        np.testing.assert_allclose(m.values[0], bin_centers())
        assert m.covered.all() and m.observed.all()

    def test_reverse_full_copy(self, golden_dir, copies):
        # Milestone B: reverse copy must produce the identical 5'->3' ramp
        track = open_bedgraph(golden_dir["bedgraph"])
        m = project_continuous(track, _copy_row(copies, 20_000), gd.CONS_LEN, BIN)
        np.testing.assert_allclose(m.values[0], bin_centers())

    def test_fragment_isolation(self, golden_dir, copies):
        # Milestone C: internal fragment [400,600) fills only bins 4-5
        track = open_bedgraph(golden_dir["bedgraph"])
        m = project_continuous(track, _copy_row(copies, 50_000), gd.CONS_LEN, BIN)
        assert m.covered[0].tolist() == [False] * 4 + [True] * 2 + [False] * 4
        np.testing.assert_allclose(m.values[0, 4:6], bin_centers()[4:6])
        assert np.isnan(m.values[0, :4]).all() and np.isnan(m.values[0, 6:]).all()

    def test_reverse_fragment(self, golden_dir, copies):
        # P3R covers consensus [700,1000) only, from a reverse genomic copy
        track = open_bedgraph(golden_dir["bedgraph"])
        m = project_continuous(track, _copy_row(copies, 40_000), gd.CONS_LEN, BIN)
        assert m.covered[0, 7:].all() and not m.covered[0, :7].any()
        np.testing.assert_allclose(m.values[0, 7:], bin_centers()[7:])

    def test_all_copies_matrix(self, golden_dir, copies):
        track = open_bedgraph(golden_dir["bedgraph"])
        m = project_continuous(track, copies, gd.CONS_LEN, BIN)
        assert m.n_copies == len(gd.COPIES)
        # per-bin element counts (all 7 copies incl. SHORT and STAR)
        expected_elements = [4, 4, 3, 3, 3, 3, 2, 3, 3, 3]
        assert m.covered.sum(axis=0).tolist() == expected_elements
        assert m.observed.sum(axis=0).tolist() == expected_elements

    def test_partial_bin_copy_value(self, golden_dir, copies):
        # SHORT covers [100,150): its bin-1 value is mean(100..149) = 124.5
        track = open_bedgraph(golden_dir["bedgraph"])
        m = project_continuous(track, _copy_row(copies, 60_000), gd.CONS_LEN, BIN)
        assert m.values[0, 1] == pytest.approx(124.5)
        assert m.covered[0].tolist() == [False, True] + [False] * 8

    def test_signal_gap_covered_but_not_observed(self, copies, tmp_path):
        # a track with data only over consensus [0,500) of copy F
        p = tmp_path / "half.bedGraph"
        p.write_text(f"{gd.CHROM}\t10000\t10500\t2.0\n")
        track = open_bedgraph(str(p))
        m = project_continuous(track, _copy_row(copies, 10_000), gd.CONS_LEN, BIN)
        assert m.covered[0].all()
        assert m.observed[0].tolist() == [True] * 5 + [False] * 5
        assert np.isnan(m.values[0, 5:]).all()
        np.testing.assert_allclose(m.values[0, :5], 2.0)


class TestIntervalProjection:
    def test_projection(self, golden_dir, copies):
        # Milestone F: peaks land at the correct consensus positions
        track = open_intervals(golden_dir["peaks_bed"])
        proj = project_intervals(track, copies)
        got = {r.peak_name: (r.cons_start, r.cons_end) for r in proj.itertuples()}
        assert got == gd.EXPECTED_PEAK_PROJECTIONS
        # peak_out overlaps nothing and must not appear
        assert "peak_out" not in set(proj["peak_name"])

    def test_occurrence(self, golden_dir, copies):
        track = open_intervals(golden_dir["peaks_bed"])
        proj = project_intervals(track, copies)
        edges = make_bins(gd.CONS_LEN, BIN)
        cov = element_coverage(copies, edges)
        occ = interval_occurrence(proj, cov, copies.index.to_numpy(), edges)
        assert occ["peak_count"].tolist() == [0, 0, 1, 0, 0, 0, 0, 0, 1, 0]
        assert occ["copies_with_peak"].tolist() == [0, 0, 1, 0, 0, 0, 0, 0, 1, 0]
        # denominators: bin 2 covered by 3 copies -> fraction 1/3
        assert occ["n_elements"].tolist() == [4, 4, 3, 3, 3, 3, 2, 3, 3, 3]
        assert occ["fraction_copies_with_peak"][2] == pytest.approx(1 / 3)
        assert occ["fraction_copies_with_peak"][8] == pytest.approx(1 / 3)

    def test_element_coverage_geometry(self, copies):
        edges = make_bins(gd.CONS_LEN, BIN)
        cov = element_coverage(copies, edges)
        assert cov.sum(axis=0).tolist() == [4, 4, 3, 3, 3, 3, 2, 3, 3, 3]


class TestPointProjection:
    def test_points(self, golden_dir, copies):
        track = open_points(golden_dir["sites"])
        pts = project_points(track, copies)
        got = dict(zip(pts["position"], pts["cons_position"]))
        assert got[10_050] == pytest.approx(50.5)
        assert got[20_950] == pytest.approx(49.5)  # reverse copy R
        assert len(pts) == 2
