"""Project continuous signal onto repeat-consensus bins, per TE copy.

For each selected copy the per-base signal over its genomic interval is
fetched, oriented 5'->3' along the consensus (strand-aware), mapped to
consensus positions (affine within the hit) and averaged into consensus
bins. The result keeps a full copies x bins matrix so aggregation can
report honest denominators (how many copies cover each bin, and how many
actually contributed data).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..chromnames import ChromReport
from ..coordinates.mapping import CopyAlignment, base_consensus_positions, make_bins
from ..tracks.base import ContinuousTrack


@dataclass
class CopyBinMatrix:
    """Per-copy binned signal on consensus coordinates for one track.

    values[i, j]  : mean signal of copy i in consensus bin j (NaN when the
                    copy does not cover the bin, or covers it with no data).
    covered[i, j] : True when copy i's consensus footprint overlaps bin j.
    observed[i, j]: True when copy i contributed >= 1 finite base value to bin j.
    """

    track_name: str
    bin_edges: np.ndarray  # length n_bins + 1
    copy_index: np.ndarray  # annotation-frame index of each row
    values: np.ndarray  # (n_copies, n_bins) float32
    covered: np.ndarray  # (n_copies, n_bins) bool
    observed: np.ndarray  # (n_copies, n_bins) bool
    n_skipped_chrom: int = 0  # copies skipped: chromosome absent from track

    @property
    def n_bins(self) -> int:
        return len(self.bin_edges) - 1

    @property
    def n_copies(self) -> int:
        return self.values.shape[0]


def copy_alignments(copies: pd.DataFrame) -> list[CopyAlignment]:
    return [
        CopyAlignment(
            genome_start=int(r.genome_start),
            genome_end=int(r.genome_end),
            cons_start=int(r.cons_start),
            cons_end=int(r.cons_end),
            strand=str(r.strand),
        )
        for r in copies.itertuples()
    ]


def project_continuous(
    track: ContinuousTrack,
    copies: pd.DataFrame,
    cons_length: int,
    bin_size: int,
    chrom_report: ChromReport | None = None,
) -> CopyBinMatrix:
    """Build the copies x bins matrix for one continuous track.

    ``copies`` is a (filtered) canonical annotation frame restricted to one
    repeat name. ``chrom_report`` maps reference chromosome names to the
    track's own naming; None means names are used as-is.
    """
    edges = make_bins(cons_length, bin_size)
    n_bins = len(edges) - 1
    n_copies = len(copies)
    values = np.full((n_copies, n_bins), np.nan, dtype="float32")
    covered = np.zeros((n_copies, n_bins), dtype=bool)
    observed = np.zeros((n_copies, n_bins), dtype=bool)
    skipped = 0

    aligns = copy_alignments(copies)
    chrom_col = copies["chrom"].astype(str).to_numpy()

    for i, aln in enumerate(aligns):
        ref_chrom = chrom_col[i]
        track_chrom = chrom_report.translate(ref_chrom) if chrom_report else ref_chrom
        if track_chrom is None:
            skipped += 1
            continue
        cons_pos = base_consensus_positions(aln)
        bin_idx = np.minimum(
            (cons_pos // bin_size).astype("int64"), n_bins - 1
        )
        bin_idx = np.clip(bin_idx, 0, n_bins - 1)
        vals = track.values(track_chrom, aln.genome_start, aln.genome_end)

        base_cov = np.bincount(bin_idx, minlength=n_bins)
        covered[i] = base_cov > 0

        finite = np.isfinite(vals)
        if finite.any():
            fin_idx = bin_idx[finite]
            sums = np.bincount(fin_idx, weights=vals[finite], minlength=n_bins)
            counts = np.bincount(fin_idx, minlength=n_bins)
            has_obs = counts > 0
            observed[i] = has_obs
            with np.errstate(invalid="ignore", divide="ignore"):
                row = sums / counts
            values[i, has_obs] = row[has_obs].astype("float32")

    return CopyBinMatrix(
        track_name=track.name,
        bin_edges=edges,
        copy_index=copies.index.to_numpy(),
        values=values,
        covered=covered,
        observed=observed,
        n_skipped_chrom=skipped,
    )
