"""Project interval (peak) and point tracks onto consensus coordinates."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..chromnames import ChromReport
from ..coordinates.mapping import (
    consensus_point_of,
    genomic_interval_to_consensus,
)
from ..tracks.base import IntervalTrack, PointTrack
from .continuous import copy_alignments

PROJECTED_INTERVAL_COLUMNS = [
    "copy_idx",  # annotation-frame index of the TE copy
    "chrom",
    "peak_start",  # original genomic peak coordinates
    "peak_end",
    "cons_start",  # projected consensus coordinates (clipped to the copy)
    "cons_end",
    "peak_name",
    "peak_score",
    "peak_strand",
]

PROJECTED_POINT_COLUMNS = [
    "copy_idx",
    "chrom",
    "position",  # genomic position of the event
    "cons_position",  # fractional consensus position
    "value",
]


def project_intervals(
    track: IntervalTrack,
    copies: pd.DataFrame,
    chrom_report: ChromReport | None = None,
) -> pd.DataFrame:
    """One row per (TE copy, overlapping peak), with projected coordinates.

    A peak overlapping several copies is reported once per copy — the
    denominator for occurrence statistics is the copy, not the peak.
    Projected coordinates are clipped to the copy's consensus footprint.
    """
    rows = []
    aligns = copy_alignments(copies)
    chrom_col = copies["chrom"].astype(str).to_numpy()
    idx = copies.index.to_numpy()
    for i, aln in enumerate(aligns):
        ref_chrom = chrom_col[i]
        track_chrom = chrom_report.translate(ref_chrom) if chrom_report else ref_chrom
        if track_chrom is None:
            continue
        for iv in track.fetch(track_chrom, aln.genome_start, aln.genome_end):
            proj = genomic_interval_to_consensus(aln, iv.start, iv.end)
            if proj is None:
                continue
            rows.append(
                (
                    idx[i],
                    ref_chrom,
                    iv.start,
                    iv.end,
                    proj[0],
                    proj[1],
                    iv.name,
                    iv.score,
                    iv.strand,
                )
            )
    return pd.DataFrame(rows, columns=PROJECTED_INTERVAL_COLUMNS)


def project_points(
    track: PointTrack,
    copies: pd.DataFrame,
    chrom_report: ChromReport | None = None,
) -> pd.DataFrame:
    """One row per (TE copy, event falling inside the copy)."""
    rows = []
    aligns = copy_alignments(copies)
    chrom_col = copies["chrom"].astype(str).to_numpy()
    idx = copies.index.to_numpy()
    for i, aln in enumerate(aligns):
        ref_chrom = chrom_col[i]
        track_chrom = chrom_report.translate(ref_chrom) if chrom_report else ref_chrom
        if track_chrom is None:
            continue
        pos, vals = track.fetch(track_chrom, aln.genome_start, aln.genome_end)
        for p, v in zip(pos, vals):
            cp = consensus_point_of(aln, int(p))
            if cp is None:
                continue
            rows.append((idx[i], ref_chrom, int(p), float(cp), float(v)))
    return pd.DataFrame(rows, columns=PROJECTED_POINT_COLUMNS)


def interval_occurrence(
    projected: pd.DataFrame,
    covered: np.ndarray,
    copy_index: np.ndarray,
    bin_edges: np.ndarray,
) -> pd.DataFrame:
    """Per-bin occurrence statistics for projected intervals.

    Two distinct quantities (different denominators, both reported):

    - ``peak_count``: number of projected (copy, peak) intervals
      overlapping the bin.
    - ``copies_with_peak`` / ``fraction_copies_with_peak``: number of TE
      copies with >= 1 peak overlapping the bin, divided by the number of
      copies whose consensus footprint covers the bin (``n_elements``).
    """
    n_bins = len(bin_edges) - 1
    starts = bin_edges[:-1]
    ends = bin_edges[1:]
    peak_count = np.zeros(n_bins, dtype="int64")
    copy_pos = {int(c): i for i, c in enumerate(copy_index)}
    copy_has_peak = np.zeros((len(copy_index), n_bins), dtype=bool)

    for row in projected.itertuples():
        b0 = int(np.searchsorted(ends, row.cons_start, side="right"))
        b1 = int(np.searchsorted(starts, row.cons_end, side="left"))
        b0 = max(0, min(b0, n_bins))
        b1 = max(0, min(b1, n_bins))
        if b1 <= b0:
            continue
        peak_count[b0:b1] += 1
        ci = copy_pos.get(int(row.copy_idx))
        if ci is not None:
            copy_has_peak[ci, b0:b1] = True

    n_elements = covered.sum(axis=0).astype("int64")
    copies_with_peak = copy_has_peak.sum(axis=0).astype("int64")
    with np.errstate(invalid="ignore", divide="ignore"):
        fraction = np.where(n_elements > 0, copies_with_peak / n_elements, np.nan)
    return pd.DataFrame(
        {
            "consensus_bin": np.arange(n_bins),
            "consensus_start": starts,
            "consensus_end": ends,
            "peak_count": peak_count,
            "copies_with_peak": copies_with_peak,
            "n_elements": n_elements,
            "fraction_copies_with_peak": fraction,
        }
    )
