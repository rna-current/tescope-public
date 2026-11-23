"""Track abstractions.

Data semantics (what a track *is*) live here; visualization semantics
(colors, panel order, visibility) are display properties carried in
:class:`TrackMeta.display` and interpreted only by the plotting/GUI layers.

```
Track
├── ContinuousTrack   per-base quantitative signal   (bigWig, bedGraph, BAM coverage)
├── IntervalTrack     intervals/peaks                (BED, narrowPeak, broadPeak)
└── PointTrack        single-base events             (site lists, BED intervals of length 1)
```
"""

from __future__ import annotations

import abc
from dataclasses import asdict, dataclass, field
from enum import Enum

import numpy as np


class TrackKind(str, Enum):
    CONTINUOUS = "continuous"
    INTERVAL = "interval"
    POINT = "point"


@dataclass
class TrackMeta:
    """User-facing metadata shared by all track kinds."""

    name: str
    path: str
    kind: TrackKind
    format: str = ""  # "bigwig" | "bedgraph" | "bed" | "narrowpeak" | "broadpeak" | "bam" | "points"
    sample: str = ""
    group: str = ""
    replicate: int | None = None
    strand: str | None = None  # optional strand the signal belongs to
    metadata: dict = field(default_factory=dict)
    display: dict = field(default_factory=dict)  # color, visible, order, ...

    def to_dict(self) -> dict:
        d = asdict(self)
        d["kind"] = self.kind.value
        return d

    @classmethod
    def from_dict(cls, d: dict) -> TrackMeta:
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        known["kind"] = TrackKind(known["kind"])
        return cls(**known)


class Track(abc.ABC):
    """A loaded, queryable signal source."""

    def __init__(self, meta: TrackMeta):
        self.meta = meta

    @property
    def name(self) -> str:
        return self.meta.name

    @abc.abstractmethod
    def chroms(self) -> dict[str, int]:
        """Chromosome name -> length (0 when the source has no lengths)."""

    def close(self) -> None:  # noqa: B027 - optional hook, not abstract
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


class ContinuousTrack(Track):
    @abc.abstractmethod
    def values(self, chrom: str, start: int, end: int) -> np.ndarray:
        """Per-base values over [start, end); NaN where no data.

        Always returns a float64 array of length ``end - start``. Queries on
        chromosomes absent from the track return all-NaN (callers decide how
        to report missingness; see chromosome-name reconciliation).
        """


@dataclass
class Interval:
    start: int
    end: int
    name: str = ""
    score: float = float("nan")
    strand: str = "."


class IntervalTrack(Track):
    @abc.abstractmethod
    def fetch(self, chrom: str, start: int, end: int) -> list[Interval]:
        """Intervals overlapping [start, end), genomic coordinates."""

    @abc.abstractmethod
    def count(self) -> int:
        """Total number of intervals in the track."""


class PointTrack(Track):
    @abc.abstractmethod
    def fetch(self, chrom: str, start: int, end: int) -> tuple[np.ndarray, np.ndarray]:
        """(positions, values) of events with position in [start, end).

        ``values`` is NaN where the source has no score.
        """

    @abc.abstractmethod
    def count(self) -> int:
        """Total number of points in the track."""
