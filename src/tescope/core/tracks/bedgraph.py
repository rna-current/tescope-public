"""bedGraph continuous tracks.

The whole file is parsed once into per-chromosome sorted interval arrays
(start, end, value); per-base queries are answered with binary search.
bedGraph files are typically far smaller than the genomes they describe,
so in-memory interval arrays are acceptable for v1 (documented in
docs/limitations.md; bigWig is recommended for genome-wide dense signal).
"""

from __future__ import annotations

import gzip
import os

import numpy as np

from .base import ContinuousTrack, TrackKind, TrackMeta


class BedGraphParseError(ValueError):
    pass


def _open_text(path: str):
    if os.fspath(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


class BedGraphTrack(ContinuousTrack):
    def __init__(self, meta: TrackMeta):
        super().__init__(meta)
        self._by_chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        self._load(meta.path)

    def _load(self, path: str) -> None:
        chroms: dict[str, list[tuple[int, int, float]]] = {}
        with _open_text(path) as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                f = line.split()
                if len(f) < 4:
                    raise BedGraphParseError(
                        f"{path}:{lineno}: expected 4 columns (chrom start end value), got {len(f)}"
                    )
                try:
                    start, end, value = int(f[1]), int(f[2]), float(f[3])
                except ValueError as e:
                    raise BedGraphParseError(f"{path}:{lineno}: {e}") from None
                if end < start:
                    raise BedGraphParseError(f"{path}:{lineno}: end < start")
                if end == start:
                    continue
                chroms.setdefault(f[0], []).append((start, end, value))
        if not chroms:
            raise BedGraphParseError(f"{path}: no data lines found")
        for chrom, rows in chroms.items():
            rows.sort()
            starts = np.array([r[0] for r in rows], dtype="int64")
            ends = np.array([r[1] for r in rows], dtype="int64")
            vals = np.array([r[2] for r in rows], dtype="float64")
            if (starts[1:] < ends[:-1]).any():
                raise BedGraphParseError(
                    f"{path}: overlapping intervals on {chrom}; "
                    "bedGraph tracks must be non-overlapping"
                )
            self._by_chrom[chrom] = (starts, ends, vals)

    def chroms(self) -> dict[str, int]:
        return {c: int(ends[-1]) for c, (_, ends, _) in self._by_chrom.items()}

    def values(self, chrom: str, start: int, end: int) -> np.ndarray:
        out = np.full(end - start, np.nan, dtype="float64")
        entry = self._by_chrom.get(chrom)
        if entry is None:
            return out
        starts, ends, vals = entry
        # candidate intervals overlapping [start, end)
        lo = int(np.searchsorted(ends, start, side="right"))
        hi = int(np.searchsorted(starts, end, side="left"))
        for i in range(lo, hi):
            a = max(int(starts[i]), start)
            b = min(int(ends[i]), end)
            if b > a:
                out[a - start : b - start] = vals[i]
        return out


def open_bedgraph(path: str, name: str | None = None, **meta_kwargs) -> BedGraphTrack:
    meta = TrackMeta(
        name=name or os.path.basename(path),
        path=os.fspath(path),
        kind=TrackKind.CONTINUOUS,
        format="bedgraph",
        **meta_kwargs,
    )
    return BedGraphTrack(meta)
