"""The per-locus table shown in the Loci panel and exported as loci.tsv.

Built from data the analysis already holds (the filtered ``copies``
frame and per-track matrices/projections). Column semantics:

- ``locus_id``: ``chrom:start-end(strand)``, 0-based half-open in the
  annotation assembly — stable across sessions for the same
  annotation input. Display widgets show 1-based coordinates; this ID
  stays 0-based because it is the machine-readable key used by
  selection files and the estimate-import contract.
- ``mean_<track>``: the copy's mean signal over its *observed* bins
  for that continuous track (input units, unchanged). NaN = the copy
  contributed no finite data (missing), which is distinct from a
  measured 0. ``obs_bins_<track>`` gives the number of observed bins
  (effective support) next to every mean.
- ``peaks_<track>``: number of projected intervals overlapping the
  copy for an interval track; 0 here is a real observed zero (the
  track covered the copy's chromosome), except when the copy's
  chromosome was absent from the track — then NaN.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..analysis import AnalysisResult
from ..tracks.base import TrackKind

#: columns copied verbatim from the filtered annotation frame
_COPY_COLUMNS = [
    "chrom",
    "genome_start",
    "genome_end",
    "strand",
    "repeat_name",
    "repeat_class",
    "repeat_family",
    "rm_id",
    "is_secondary",
    "genome_length",
    "cons_start",
    "cons_end",
    "consensus_coverage",
    "completeness",
]


def locus_id_for_row(row) -> str:
    """Stable locus ID for one annotation row (0-based half-open)."""
    return (
        f"{row.chrom}:{int(row.genome_start)}-{int(row.genome_end)}"
        f"({row.strand})"
    )


def locus_ids(copies: pd.DataFrame) -> pd.Series:
    return (
        copies["chrom"].astype(str)
        + ":"
        + copies["genome_start"].astype(int).astype(str)
        + "-"
        + copies["genome_end"].astype(int).astype(str)
        + "("
        + copies["strand"].astype(str)
        + ")"
    )


def build_locus_table(result: AnalysisResult) -> pd.DataFrame:
    """One row per analyzed copy, indexed like ``result.copies``."""
    copies = result.copies
    cols = [c for c in _COPY_COLUMNS if c in copies.columns]
    table = copies[cols].copy()
    table.insert(0, "locus_id", locus_ids(copies))

    # rm_id fragment linkage: how many analyzed records share this rm_id
    if "rm_id" in table.columns:
        frag_counts = table.groupby("rm_id")["rm_id"].transform("size")
        table["rm_id_fragments"] = frag_counts.astype("int64")

    for t in result.tracks:
        label = t.name
        if t.kind is TrackKind.CONTINUOUS and t.matrix is not None:
            m = t.matrix
            vals = m.values.astype("float64")
            obs = np.isfinite(vals)
            n_obs = obs.sum(axis=1).astype("int64")
            with np.errstate(invalid="ignore"):
                mean = np.where(
                    n_obs > 0,
                    np.nansum(np.where(obs, vals, 0.0), axis=1)
                    / np.maximum(n_obs, 1),
                    np.nan,
                )
            per_copy = pd.Series(mean, index=pd.Index(m.copy_index))
            per_obs = pd.Series(n_obs, index=pd.Index(m.copy_index))
            table[f"mean_{label}"] = per_copy.reindex(table.index)
            table[f"obs_bins_{label}"] = (
                per_obs.reindex(table.index).fillna(0).astype("int64")
            )
        elif t.kind is TrackKind.INTERVAL and t.projected is not None:
            counts = (
                t.projected.groupby("copy_idx").size()
                if len(t.projected)
                else pd.Series(dtype="int64")
            )
            col = counts.reindex(table.index).fillna(0).astype("float64")
            # chromosome-absent copies: evidence unavailable, not zero.
            # The matrix records this for continuous tracks; for interval
            # tracks we mark copies whose chromosome the track could not
            # translate (conservative: only when the report says so).
            if t.chrom_report is not None:
                report = t.chrom_report
                untranslatable = ~copies["chrom"].astype(str).map(
                    lambda c, _r=report: _r.translate(c) is not None
                )
                col[untranslatable.to_numpy()] = np.nan
            table[f"peaks_{label}"] = col
    return table


def signal_columns(table: pd.DataFrame) -> list[str]:
    return [
        c
        for c in table.columns
        if c.startswith(("mean_", "peaks_", "est_"))
    ]
