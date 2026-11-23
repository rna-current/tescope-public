"""BED / narrowPeak / broadPeak interval tracks and BED-derived point tracks.

Intervals are loaded into per-chromosome sorted arrays and queried by
binary search. Peak files remain interval/event tracks; they are never
converted into quantitative curves.
"""

from __future__ import annotations

import gzip
import os

import numpy as np

from .base import Interval, IntervalTrack, PointTrack, TrackKind, TrackMeta

_FORMAT_MIN_COLS = {"bed": 3, "narrowpeak": 10, "broadpeak": 9}


class BedParseError(ValueError):
    pass


def _open_text(path: str):
    if os.fspath(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def detect_bed_format(path: str) -> str:
    p = os.fspath(path).lower()
    for ext, fmt in (
        (".narrowpeak", "narrowpeak"),
        (".broadpeak", "broadpeak"),
        (".bed", "bed"),
    ):
        if p.endswith(ext) or p.endswith(ext + ".gz"):
            return fmt
    return "bed"


def _parse_bed_like(path: str, fmt: str):
    """Returns dict chrom -> (starts, ends, names, scores, strands)."""
    min_cols = _FORMAT_MIN_COLS[fmt]
    per_chrom: dict[str, list[tuple]] = {}
    with _open_text(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) == 1:  # tolerate space-separated files
                f = line.split()
            if len(f) < min_cols:
                raise BedParseError(
                    f"{path}:{lineno}: {fmt} requires >= {min_cols} columns, got {len(f)}"
                )
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError as e:
                raise BedParseError(f"{path}:{lineno}: {e}") from None
            if end <= start:
                raise BedParseError(f"{path}:{lineno}: end <= start")
            name = f[3] if len(f) > 3 else ""
            score = np.nan
            if len(f) > 4 and f[4] not in (".", ""):
                try:
                    score = float(f[4])
                except ValueError:
                    score = np.nan
            strand = f[5] if len(f) > 5 and f[5] in ("+", "-") else "."
            # narrowPeak signalValue (col 7) is a better magnitude than the
            # capped BED score; expose it as the score for narrow/broad peaks.
            if fmt in ("narrowpeak", "broadpeak") and len(f) > 6:
                try:
                    score = float(f[6])
                except ValueError:
                    pass
            per_chrom.setdefault(f[0], []).append((start, end, name, score, strand))
    if not per_chrom:
        raise BedParseError(f"{path}: no intervals found")
    out = {}
    for chrom, rows in per_chrom.items():
        rows.sort(key=lambda r: (r[0], r[1]))
        out[chrom] = (
            np.array([r[0] for r in rows], dtype="int64"),
            np.array([r[1] for r in rows], dtype="int64"),
            [r[2] for r in rows],
            np.array([r[3] for r in rows], dtype="float64"),
            [r[4] for r in rows],
        )
    return out


class BedIntervalTrack(IntervalTrack):
    def __init__(self, meta: TrackMeta):
        super().__init__(meta)
        fmt = meta.format or detect_bed_format(meta.path)
        self._data = _parse_bed_like(meta.path, fmt)
        self._max_len = {
            c: int((ends - starts).max()) for c, (starts, ends, *_rest) in self._data.items()
        }
        self._n = sum(len(v[0]) for v in self._data.values())

    def chroms(self) -> dict[str, int]:
        return {c: int(v[1][-1]) for c, v in self._data.items()}

    def count(self) -> int:
        return self._n

    def fetch(self, chrom: str, start: int, end: int) -> list[Interval]:
        entry = self._data.get(chrom)
        if entry is None:
            return []
        starts, ends, names, scores, strands = entry
        # sorted by start; an overlapping interval starts within
        # [start - max_len, end)
        lo = int(np.searchsorted(starts, start - self._max_len[chrom], side="left"))
        hi = int(np.searchsorted(starts, end, side="left"))
        out = []
        for i in range(lo, hi):
            if ends[i] > start:
                out.append(
                    Interval(
                        start=int(starts[i]),
                        end=int(ends[i]),
                        name=names[i],
                        score=float(scores[i]),
                        strand=strands[i],
                    )
                )
        return out


class BedPointTrack(PointTrack):
    """Point events from a BED-like file.

    Each interval contributes its single base when length == 1, otherwise
    its midpoint (reported once, at ``(start + end) // 2``).
    """

    def __init__(self, meta: TrackMeta):
        super().__init__(meta)
        data = _parse_bed_like(meta.path, "bed")
        self._points: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        self._n = 0
        for chrom, (starts, ends, _names, scores, _strands) in data.items():
            pos = np.where(ends - starts == 1, starts, (starts + ends) // 2)
            order = np.argsort(pos, kind="stable")
            self._points[chrom] = (pos[order], scores[order])
            self._n += len(pos)

    def chroms(self) -> dict[str, int]:
        return {c: int(pos[-1]) + 1 for c, (pos, _) in self._points.items()}

    def count(self) -> int:
        return self._n

    def fetch(self, chrom: str, start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
        entry = self._points.get(chrom)
        if entry is None:
            return np.empty(0, dtype="int64"), np.empty(0, dtype="float64")
        pos, vals = entry
        lo = int(np.searchsorted(pos, start, side="left"))
        hi = int(np.searchsorted(pos, end, side="left"))
        return pos[lo:hi], vals[lo:hi]


def open_intervals(path: str, name: str | None = None, **meta_kwargs) -> BedIntervalTrack:
    fmt = meta_kwargs.pop("format", "") or detect_bed_format(path)
    meta = TrackMeta(
        name=name or os.path.basename(path),
        path=os.fspath(path),
        kind=TrackKind.INTERVAL,
        format=fmt,
        **meta_kwargs,
    )
    return BedIntervalTrack(meta)


def open_points(path: str, name: str | None = None, **meta_kwargs) -> BedPointTrack:
    meta = TrackMeta(
        name=name or os.path.basename(path),
        path=os.fspath(path),
        kind=TrackKind.POINT,
        format="points",
        **meta_kwargs,
    )
    return BedPointTrack(meta)
