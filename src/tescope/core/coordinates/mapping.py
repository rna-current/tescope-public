"""Genomic <-> repeat-consensus coordinate mapping.

Each RepeatMasker hit aligns a genomic interval ``[genome_start, genome_end)``
to a consensus interval ``[cons_start, cons_end)``. The ``.out`` format does
not retain the per-base alignment, so TEscope maps positions **linearly
(affine) within each hit**: insertions/deletions inside a copy are
distributed proportionally. This is the standard approximation for
RepeatMasker-based consensus projections and is documented as a limitation.

Orientation: consensus coordinates always run 5'->3' along the repeat
model. For a hit on the genomic ``-`` strand, the *last* genomic base of
the interval corresponds to the *first* consensus base of the hit
(``cons_start``), and signal arrays must be reversed before binning.

All functions treat coordinates as 0-based half-open.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CopyAlignment:
    """Minimal alignment geometry of one TE copy (one RepeatMasker hit)."""

    genome_start: int
    genome_end: int
    cons_start: int
    cons_end: int
    strand: str  # '+' or '-'

    def __post_init__(self):
        if self.genome_end <= self.genome_start:
            raise ValueError("empty genomic interval")
        if self.cons_end < self.cons_start:
            raise ValueError("cons_end < cons_start")
        if self.strand not in ("+", "-"):
            raise ValueError(f"invalid strand {self.strand!r}")

    @property
    def genome_length(self) -> int:
        return self.genome_end - self.genome_start

    @property
    def cons_span(self) -> int:
        return self.cons_end - self.cons_start

    @property
    def scale(self) -> float:
        """Consensus bases per genomic base within this hit."""
        return self.cons_span / self.genome_length


def base_consensus_positions(copy: CopyAlignment) -> np.ndarray:
    """Fractional consensus position of every genomic base of the copy.

    Returns an array of length ``genome_length``; element ``i`` is the
    consensus position (float) of genomic base ``genome_start + i``. Base
    centers are used so that a hit whose consensus span differs from its
    genomic length (indels) distributes bases uniformly without edge bias.

    For '-' strand copies the array still corresponds to genomic order
    (index 0 = leftmost genomic base); values decrease along the array.
    """
    L = copy.genome_length
    offsets = np.arange(L, dtype="float64") + 0.5
    if copy.strand == "-":
        offsets = L - offsets  # reverse: last genomic base -> consensus start
    return copy.cons_start + offsets * copy.scale


def genomic_interval_to_consensus(
    copy: CopyAlignment, start: int, end: int
) -> tuple[int, int] | None:
    """Project a genomic interval onto the copy's consensus coordinates.

    The interval is first clipped to the copy; ``None`` is returned when
    there is no overlap. The result is an integer 0-based half-open
    consensus interval of length >= 1.
    """
    a = max(start, copy.genome_start)
    b = min(end, copy.genome_end)
    if b <= a:
        return None
    s = copy.scale
    if copy.strand == "+":
        c0 = copy.cons_start + (a - copy.genome_start) * s
        c1 = copy.cons_start + (b - copy.genome_start) * s
    else:
        c0 = copy.cons_start + (copy.genome_end - b) * s
        c1 = copy.cons_start + (copy.genome_end - a) * s
    lo = int(np.floor(c0))
    hi = int(np.ceil(c1))
    if hi <= lo:
        hi = lo + 1
    # clip to the copy's consensus footprint
    lo = max(lo, copy.cons_start)
    hi = min(hi, max(copy.cons_end, lo + 1))
    return lo, hi


def consensus_interval_to_genomic(
    copy: CopyAlignment, cons_start: int, cons_end: int
) -> tuple[int, int, bool] | None:
    """Project a consensus interval onto a copy's genomic coordinates.

    Exact inverse of :func:`genomic_interval_to_consensus` under the
    same per-hit affine approximation — the ``.out`` format retains no
    per-base alignment, so this is NOT gap-aware or base-accurate and
    callers must label the result "projected". Returns
    ``(genome_start, genome_end, partial)`` 0-based half-open, where
    ``partial`` is True when the feature extends beyond the copy's
    consensus footprint (only the overlapping part is projected), or
    ``None`` when the feature lies entirely outside the footprint
    (no verified correspondence — nothing may be drawn).
    """
    a = max(cons_start, copy.cons_start)
    b = min(cons_end, copy.cons_end)
    if b <= a:
        return None
    partial = a > cons_start or b < cons_end
    inv = 1.0 / copy.scale if copy.scale else 0.0
    if copy.strand == "+":
        g0 = copy.genome_start + (a - copy.cons_start) * inv
        g1 = copy.genome_start + (b - copy.cons_start) * inv
    else:
        g0 = copy.genome_end - (b - copy.cons_start) * inv
        g1 = copy.genome_end - (a - copy.cons_start) * inv
    lo = int(np.floor(g0))
    hi = int(np.ceil(g1))
    if hi <= lo:
        hi = lo + 1
    lo = max(lo, copy.genome_start)
    hi = min(hi, max(copy.genome_end, lo + 1))
    return lo, hi, partial


def consensus_point_of(copy: CopyAlignment, pos: int) -> float | None:
    """Fractional consensus position of a single genomic base (its center)."""
    if not (copy.genome_start <= pos < copy.genome_end):
        return None
    if copy.strand == "+":
        off = pos - copy.genome_start + 0.5
    else:
        off = copy.genome_end - pos - 0.5
    return copy.cons_start + off * copy.scale


def make_bins(cons_length: int, bin_size: int) -> np.ndarray:
    """Bin edges over ``[0, cons_length)`` with the last bin possibly short.

    Returns an array of edges of length ``n_bins + 1``.
    """
    if cons_length <= 0:
        raise ValueError("consensus length must be positive")
    if bin_size <= 0:
        raise ValueError("bin size must be positive")
    edges = np.arange(0, cons_length, bin_size, dtype="int64")
    return np.append(edges, cons_length)
