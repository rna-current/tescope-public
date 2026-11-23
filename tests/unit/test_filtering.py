"""TE filtering: length, consensus coverage, completeness classification."""

from __future__ import annotations

import golden_dataset as gd
import pytest

from tescope.core.annotation.repeatmasker import parse_repeatmasker_out
from tescope.core.annotation.schema import add_derived_columns
from tescope.core.filtering.filters import (
    FilterConfig,
    apply_filters,
    classify_completeness,
)


@pytest.fixture(scope="module")
def tst1(golden_dir):
    df = parse_repeatmasker_out(golden_dir["out"])
    df = df[df["repeat_name"] == "TST1"]
    return add_derived_columns(df)


def _by_start(df):
    return {int(r.genome_start): r for r in df.itertuples()}


class TestDerived:
    def test_lengths_and_coverage(self, tst1):
        rows = _by_start(tst1)
        assert rows[10_000].genome_length == 1000
        assert rows[10_000].consensus_coverage == pytest.approx(1.0)
        assert rows[30_000].consensus_coverage == pytest.approx(0.4)
        assert rows[40_000].consensus_coverage == pytest.approx(0.3)
        assert rows[60_000].consensus_coverage == pytest.approx(0.05)


class TestCompleteness:
    def test_expected_classes(self, tst1):
        cls = classify_completeness(tst1)
        got = dict(zip(tst1["genome_start"], cls.astype(str)))
        for name, (gs, *_rest) in gd.COPIES.items():
            assert got[gs] == gd.EXPECTED_COMPLETENESS[name], name

    def test_tolerance_matters(self, tst1):
        # with a huge tolerance everything reaching either half becomes full
        cls = classify_completeness(tst1, end_tolerance_frac=0.7)
        assert (cls.astype(str) == "full_length").sum() >= 4


class TestFilters:
    def test_no_filters_keeps_all(self, tst1):
        out, report = apply_filters(tst1, FilterConfig())
        assert len(out) == len(tst1)
        assert report.n_input == len(tst1)
        assert report.n_remaining == len(tst1)

    def test_min_length_excludes_short(self, tst1):
        out, report = apply_filters(tst1, FilterConfig(min_genome_length=60))
        assert 60_000 not in set(out["genome_start"])
        assert len(out) == len(tst1) - 1
        assert report.steps[-1].n_remaining == len(out)

    def test_consensus_coverage(self, tst1):
        out, _ = apply_filters(tst1, FilterConfig(min_consensus_coverage=0.5))
        assert set(out["genome_start"]) == {10_000, 20_000}

    def test_completeness_filter(self, tst1):
        out, _ = apply_filters(tst1, FilterConfig(completeness=["full_length"]))
        assert set(out["genome_start"]) == {10_000, 20_000}
        out, _ = apply_filters(
            tst1, FilterConfig(completeness=["five_prime_truncated", "internal_fragment"])
        )
        assert set(out["genome_start"]) == {40_000, 50_000, 60_000}

    def test_secondary_exclusion(self, tst1):
        out, _ = apply_filters(tst1, FilterConfig(include_secondary=False))
        assert 90_000 not in set(out["genome_start"])

    def test_strand_filter(self, tst1):
        out, _ = apply_filters(tst1, FilterConfig(strands=["-"]))
        assert set(out["genome_start"]) == {20_000, 40_000}

    def test_sequential_report(self, tst1):
        _, report = apply_filters(
            tst1, FilterConfig(min_genome_length=60, min_consensus_coverage=0.5)
        )
        labels = [s.label for s in report.steps]
        assert labels == ["Length >= 60 bp", "Consensus coverage >= 50%"]
        assert [s.n_remaining for s in report.steps] == [6, 2]

    def test_unknown_completeness_rejected(self):
        with pytest.raises(ValueError, match="completeness"):
            FilterConfig.from_dict({"completeness": ["fullish"]})
