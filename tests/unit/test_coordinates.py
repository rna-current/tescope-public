"""Consensus coordinate mapping and strand reversal (spec §12: dedicated tests)."""

from __future__ import annotations

import numpy as np
import pytest

from tescope.core.coordinates.mapping import (
    CopyAlignment,
    base_consensus_positions,
    consensus_point_of,
    genomic_interval_to_consensus,
    make_bins,
)

FWD = CopyAlignment(genome_start=10_000, genome_end=11_000, cons_start=0, cons_end=1000, strand="+")
REV = CopyAlignment(genome_start=20_000, genome_end=21_000, cons_start=0, cons_end=1000, strand="-")
FRAG_REV = CopyAlignment(
    genome_start=40_000, genome_end=40_300, cons_start=700, cons_end=1000, strand="-"
)


class TestBasePositions:
    def test_forward_identity(self):
        pos = base_consensus_positions(FWD)
        assert pos[0] == pytest.approx(0.5)
        assert pos[-1] == pytest.approx(999.5)
        assert np.all(np.diff(pos) > 0)

    def test_reverse_orientation(self):
        # last genomic base -> consensus start
        pos = base_consensus_positions(REV)
        assert pos[0] == pytest.approx(999.5)
        assert pos[-1] == pytest.approx(0.5)
        assert np.all(np.diff(pos) < 0)

    def test_reverse_fragment_offsets(self):
        pos = base_consensus_positions(FRAG_REV)
        assert pos[-1] == pytest.approx(700.5)
        assert pos[0] == pytest.approx(999.5)

    def test_indel_scaling(self):
        # 100 genomic bases aligned to 50 consensus bases: scale 0.5
        aln = CopyAlignment(genome_start=0, genome_end=100, cons_start=100, cons_end=150, strand="+")
        pos = base_consensus_positions(aln)
        assert len(pos) == 100
        assert pos[0] == pytest.approx(100.25)
        assert pos[-1] == pytest.approx(149.75)


class TestIntervalProjection:
    def test_forward(self):
        assert genomic_interval_to_consensus(FWD, 10_250, 10_300) == (250, 300)

    def test_reverse(self):
        assert genomic_interval_to_consensus(REV, 20_100, 20_200) == (800, 900)

    def test_clipping_to_copy(self):
        assert genomic_interval_to_consensus(FWD, 9_000, 10_100) == (0, 100)

    def test_no_overlap(self):
        assert genomic_interval_to_consensus(FWD, 70_000, 70_100) is None

    def test_minimum_width(self):
        # sub-base overlap still yields a 1 bp consensus interval
        aln = CopyAlignment(genome_start=0, genome_end=1000, cons_start=0, cons_end=10, strand="+")
        lo, hi = genomic_interval_to_consensus(aln, 500, 501)
        assert hi - lo >= 1

    def test_reverse_fragment(self):
        # genomic left edge of a reverse fragment = consensus 3' end
        assert genomic_interval_to_consensus(FRAG_REV, 40_000, 40_050) == (950, 1000)


class TestPointMapping:
    def test_forward(self):
        assert consensus_point_of(FWD, 10_050) == pytest.approx(50.5)

    def test_reverse(self):
        assert consensus_point_of(REV, 20_950) == pytest.approx(49.5)

    def test_outside(self):
        assert consensus_point_of(FWD, 11_000) is None


class TestBins:
    def test_even(self):
        edges = make_bins(1000, 100)
        assert len(edges) == 11
        assert edges[0] == 0 and edges[-1] == 1000

    def test_ragged_last_bin(self):
        edges = make_bins(1050, 100)
        assert edges[-2:].tolist() == [1000, 1050]

    def test_invalid(self):
        with pytest.raises(ValueError):
            make_bins(0, 100)
        with pytest.raises(ValueError):
            make_bins(100, 0)


class TestAlignmentValidation:
    def test_rejects_bad_strand(self):
        with pytest.raises(ValueError):
            CopyAlignment(0, 10, 0, 10, "C")

    def test_rejects_empty_interval(self):
        with pytest.raises(ValueError):
            CopyAlignment(10, 10, 0, 10, "+")
