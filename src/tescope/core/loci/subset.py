"""Explicit subset analysis over an existing result.

Recomputes per-bin aggregates for a user-selected locus subset FROM
THE ORIGINAL PER-COPY VALUES the analysis already holds (the
copies x bins matrices and projected intervals). Nothing is
re-projected, re-estimated, or reassigned; the eligible population's
aggregates are untouched and remain the comparison context.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..aggregation.aggregate import aggregate_profile
from ..analysis import AnalysisResult
from ..projection.continuous import CopyBinMatrix
from ..projection.intervals import interval_occurrence
from ..tracks.base import TrackKind
from .table import locus_ids


class SubsetError(ValueError):
    """User-facing subset failure (e.g. stale selection)."""


@dataclass
class SubsetResult:
    """Aggregates for an explicitly applied locus subset."""

    locus_ids: list[str]
    n_subset: int
    n_eligible: int
    profiles: dict[str, pd.DataFrame] = field(default_factory=dict)
    occurrences: dict[str, pd.DataFrame] = field(default_factory=dict)
    invalid_ids: list[str] = field(default_factory=list)

    @property
    def label(self) -> str:
        return f"subset (n={self.n_subset} of {self.n_eligible})"


def _restrict_matrix(m: CopyBinMatrix, mask: np.ndarray) -> CopyBinMatrix:
    return CopyBinMatrix(
        track_name=m.track_name,
        bin_edges=m.bin_edges,
        copy_index=m.copy_index[mask],
        values=m.values[mask],
        covered=m.covered[mask],
        observed=m.observed[mask],
        n_skipped_chrom=0,
    )


def subset_result(
    result: AnalysisResult, selected: list[str]
) -> SubsetResult:
    """Aggregate the selected loci. ``selected`` = stable locus IDs.

    IDs that no longer exist in the eligible population are returned
    in ``invalid_ids`` (reported, never silently dropped); at least
    two valid loci are required for an aggregate.
    """
    ids = locus_ids(result.copies)
    id_to_pos = {s: i for i, s in enumerate(ids)}
    valid, invalid = [], []
    for s in selected:
        (valid if s in id_to_pos else invalid).append(s)
    if len(valid) < 2:
        raise SubsetError(
            f"subset needs >= 2 loci in the eligible population "
            f"({len(valid)} valid, {len(invalid)} unknown/filtered-out)"
        )
    pos = np.array([id_to_pos[s] for s in valid], dtype="int64")
    mask = np.zeros(len(result.copies), dtype=bool)
    mask[pos] = True
    sub_index = result.copies.index.to_numpy()[mask]

    out = SubsetResult(
        locus_ids=valid,
        n_subset=len(valid),
        n_eligible=len(result.copies),
        invalid_ids=invalid,
    )
    for t in result.tracks:
        if t.kind is TrackKind.CONTINUOUS and t.matrix is not None:
            # matrix rows follow result.copies order (same construction)
            out.profiles[t.name] = aggregate_profile(
                _restrict_matrix(t.matrix, mask)
            )
        elif t.kind is TrackKind.INTERVAL and t.projected is not None:
            sub_proj = t.projected[t.projected["copy_idx"].isin(sub_index)]
            out.occurrences[t.name] = interval_occurrence(
                sub_proj,
                result.element_covered[mask],
                sub_index,
                result.bin_edges,
            )
    return out
