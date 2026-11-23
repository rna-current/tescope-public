"""Aggregation statistics, denominators, replicate combining, normalization."""

from __future__ import annotations

import golden_dataset as gd
import numpy as np
import pandas as pd
import pytest

from tescope.core.aggregation.aggregate import (
    aggregate_profile,
    combine_replicates,
    difference_profile,
    normalize_values,
)
from tescope.core.annotation.repeatmasker import parse_repeatmasker_out
from tescope.core.annotation.schema import add_derived_columns
from tescope.core.projection.continuous import project_continuous
from tescope.core.tracks.bedgraph import open_bedgraph

BIN = 100


@pytest.fixture(scope="module")
def profile_and_matrix(golden_dir):
    df = parse_repeatmasker_out(golden_dir["out"])
    df = add_derived_columns(df[df["repeat_name"] == "TST1"])
    track = open_bedgraph(golden_dir["bedgraph"])
    m = project_continuous(track, df, gd.CONS_LEN, BIN)
    return aggregate_profile(m), m


class TestAggregateProfile:
    def test_denominators(self, profile_and_matrix):
        prof, _ = profile_and_matrix
        assert prof["n_elements"].tolist() == [4, 4, 3, 3, 3, 3, 2, 3, 3, 3]
        assert prof["n_observations"].tolist() == [4, 4, 3, 3, 3, 3, 2, 3, 3, 3]

    def test_mean_values(self, profile_and_matrix):
        # every full-bin copy contributes the bin center; SHORT's partial bin 1
        # contributes mean(100..149)=124.5 -> bin 1 mean = (3*149.5+124.5)/4
        prof, _ = profile_and_matrix
        centers = np.arange(10) * BIN + (BIN - 1) / 2
        expected = centers.copy()
        expected[1] = (3 * 149.5 + 124.5) / 4
        np.testing.assert_allclose(prof["mean"], expected)

    def test_median_values(self, profile_and_matrix):
        prof, _ = profile_and_matrix
        centers = np.arange(10) * BIN + (BIN - 1) / 2
        np.testing.assert_allclose(prof["median"], centers)  # median resists SHORT

    def test_columns(self, profile_and_matrix):
        prof, _ = profile_and_matrix
        assert list(prof.columns) == [
            "consensus_bin", "consensus_start", "consensus_end",
            "mean", "median", "sem", "n_elements", "n_observations",
        ]

    def test_empty_bins_are_nan_with_zero_n(self, golden_dir):
        df = parse_repeatmasker_out(golden_dir["out"])
        df = add_derived_columns(df[df["repeat_name"] == "TST1"])
        only_mid = df[df["genome_start"] == 50_000]
        track = open_bedgraph(golden_dir["bedgraph"])
        m = project_continuous(track, only_mid, gd.CONS_LEN, BIN)
        prof = aggregate_profile(m)
        assert prof["n_elements"].tolist() == [0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
        assert prof["mean"].isna().tolist() == [True] * 4 + [False] * 2 + [True] * 4


class TestReplicates:
    def _fake_profile(self, values):
        n = len(values)
        return pd.DataFrame(
            {
                "consensus_bin": np.arange(n),
                "consensus_start": np.arange(n) * BIN,
                "consensus_end": (np.arange(n) + 1) * BIN,
                "mean": np.asarray(values, dtype=float),
                "median": np.asarray(values, dtype=float),
                "sem": np.zeros(n),
                "n_elements": np.full(n, 5),
                "n_observations": np.full(n, 5),
            }
        )

    def test_combine_mean_and_sem(self):
        g = combine_replicates([self._fake_profile([1, 2]), self._fake_profile([3, 4])])
        np.testing.assert_allclose(g["mean"], [2.0, 3.0])
        np.testing.assert_allclose(g["sem"], [1.0, 1.0])  # std([1,3])/sqrt(2) = 1
        assert g["n_replicates"].tolist() == [2, 2]

    def test_single_replicate_sem_nan(self):
        g = combine_replicates([self._fake_profile([1, 2])])
        assert g["sem"].isna().all()

    def test_difference(self):
        a = self._fake_profile([5, 5])
        b = self._fake_profile([2, 3])
        d = difference_profile(a, b)
        np.testing.assert_allclose(d["difference"], [3.0, 2.0])


class TestNormalization:
    def test_zscore(self):
        v = normalize_values(np.array([1.0, 2.0, 3.0]), "zscore")
        assert v[1] == pytest.approx(0.0)
        assert np.std(v, ddof=1) == pytest.approx(1.0)

    def test_minmax(self):
        v = normalize_values(np.array([2.0, 4.0, 6.0]), "minmax")
        np.testing.assert_allclose(v, [0.0, 0.5, 1.0])

    def test_none_identity(self):
        v = np.array([1.0, np.nan, 3.0])
        np.testing.assert_allclose(normalize_values(v, "none"), v, equal_nan=True)

    def test_nan_preserved(self):
        v = normalize_values(np.array([1.0, np.nan, 3.0]), "minmax")
        assert np.isnan(v[1])
        np.testing.assert_allclose(v[[0, 2]], [0.0, 1.0])

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            normalize_values(np.array([1.0]), "quantile")
