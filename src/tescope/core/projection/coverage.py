"""Geometry-only element coverage of consensus bins (no signal involved)."""

from __future__ import annotations

import numpy as np
import pandas as pd


def element_coverage(copies: pd.DataFrame, bin_edges: np.ndarray) -> np.ndarray:
    """Boolean (n_copies, n_bins) matrix: copy footprint overlaps bin.

    A copy covers bin j when its consensus interval [cons_start, cons_end)
    overlaps [edges[j], edges[j+1]). Zero-span hits (cons_start == cons_end)
    cover the single bin containing cons_start.
    """
    n_bins = len(bin_edges) - 1
    starts = bin_edges[:-1]
    ends = bin_edges[1:]
    covered = np.zeros((len(copies), n_bins), dtype=bool)
    cs = copies["cons_start"].to_numpy()
    ce = copies["cons_end"].to_numpy()
    for i in range(len(copies)):
        b0 = int(np.searchsorted(ends, cs[i], side="right"))
        b1 = int(np.searchsorted(starts, max(ce[i], cs[i] + 1), side="left"))
        b0 = max(0, min(b0, n_bins))
        b1 = max(b0, min(b1, n_bins))
        covered[i, b0:b1] = True
    return covered
