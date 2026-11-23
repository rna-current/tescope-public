"""Aggregate per-copy binned signal across TE copies, with honest denominators.

Fragmented copies mean each consensus bin is supported by a different
number of elements; every aggregate row therefore carries ``n_elements``
(copies whose footprint covers the bin) and ``n_observations`` (copies that
contributed finite data). Normalizations are applied only when explicitly
requested and are recorded by the caller in the project config.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from ..projection.continuous import CopyBinMatrix

PROFILE_COLUMNS = [
    "consensus_bin",
    "consensus_start",
    "consensus_end",
    "mean",
    "median",
    "sem",
    "n_elements",
    "n_observations",
]

VALID_STATISTICS = ("mean", "median")
VALID_NORMALIZATIONS = ("none", "zscore", "minmax")


def aggregate_profile(matrix: CopyBinMatrix) -> pd.DataFrame:
    """Aggregate a copies x bins matrix into one per-bin profile table."""
    vals = matrix.values.astype("float64")
    n_obs = np.isfinite(vals).sum(axis=0).astype("int64")
    n_el = matrix.covered.sum(axis=0).astype("int64")
    # bins with 0-1 observations legitimately produce NaN mean/sem; the numpy
    # empty-slice / ddof warnings are expected and suppressed here only
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(np.where(np.isfinite(vals), vals, np.nan), axis=0)
        median = np.nanmedian(vals, axis=0)
        std = np.nanstd(vals, axis=0, ddof=1)
    mean[n_obs == 0] = np.nan
    median[n_obs == 0] = np.nan
    sem = np.where(n_obs > 1, std / np.sqrt(np.maximum(n_obs, 1)), np.nan)
    return pd.DataFrame(
        {
            "consensus_bin": np.arange(matrix.n_bins),
            "consensus_start": matrix.bin_edges[:-1],
            "consensus_end": matrix.bin_edges[1:],
            "mean": mean,
            "median": median,
            "sem": sem,
            "n_elements": n_el,
            "n_observations": n_obs,
        }
    )


def combine_replicates(profiles: list[pd.DataFrame], statistic: str = "mean") -> pd.DataFrame:
    """Combine replicate profiles of one group into a group profile.

    The group value per bin is the mean (or median) of the replicates'
    ``statistic`` values; ``sem`` is the between-replicate SEM. Denominators
    are summarized conservatively: ``n_elements``/``n_observations`` report
    the minimum across replicates (identical for same-assembly runs).
    """
    if statistic not in VALID_STATISTICS:
        raise ValueError(f"statistic must be one of {VALID_STATISTICS}")
    if not profiles:
        raise ValueError("no profiles to combine")
    base = profiles[0]
    stacked = np.stack([p[statistic].to_numpy(dtype="float64") for p in profiles])
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        value = np.nanmean(stacked, axis=0)
        n_rep = np.isfinite(stacked).sum(axis=0)
        std = np.nanstd(stacked, axis=0, ddof=1)
    sem = np.where(n_rep > 1, std / np.sqrt(np.maximum(n_rep, 1)), np.nan)
    out = pd.DataFrame(
        {
            "consensus_bin": base["consensus_bin"],
            "consensus_start": base["consensus_start"],
            "consensus_end": base["consensus_end"],
            statistic: value,
            "sem": sem,
            "n_replicates": n_rep,
            "n_elements": np.minimum.reduce([p["n_elements"].to_numpy() for p in profiles]),
            "n_observations": np.minimum.reduce(
                [p["n_observations"].to_numpy() for p in profiles]
            ),
        }
    )
    return out


def difference_profile(a: pd.DataFrame, b: pd.DataFrame, statistic: str = "mean") -> pd.DataFrame:
    """Per-bin difference (a - b) of two group/track profiles."""
    if len(a) != len(b):
        raise ValueError("profiles have different numbers of bins")
    out = a[["consensus_bin", "consensus_start", "consensus_end"]].copy()
    out["difference"] = a[statistic].to_numpy() - b[statistic].to_numpy()
    out["n_elements"] = np.minimum(a["n_elements"].to_numpy(), b["n_elements"].to_numpy())
    return out


def normalize_values(values: np.ndarray, method: str) -> np.ndarray:
    """Explicit per-profile normalization for display ('zscore' | 'minmax' | 'none')."""
    if method not in VALID_NORMALIZATIONS:
        raise ValueError(f"normalization must be one of {VALID_NORMALIZATIONS}")
    v = np.asarray(values, dtype="float64")
    if method == "none":
        return v
    finite = np.isfinite(v)
    if finite.sum() < 2:
        return v
    if method == "zscore":
        mu = np.nanmean(v)
        sd = np.nanstd(v, ddof=1)
        return (v - mu) / sd if sd > 0 else v - mu
    lo = np.nanmin(v)
    hi = np.nanmax(v)
    return (v - lo) / (hi - lo) if hi > lo else np.where(finite, 0.0, np.nan)
