"""Open a track file by sniffing its format from extension/content."""

from __future__ import annotations

import os

from .base import Track, TrackKind
from .bed import detect_bed_format, open_intervals, open_points
from .bedgraph import open_bedgraph
from .bigwig import open_bigwig

_CONTINUOUS_EXT = (".bw", ".bigwig")
_BEDGRAPH_EXT = (".bedgraph", ".bg", ".bdg")
_INTERVAL_EXT = (".bed", ".narrowpeak", ".broadpeak")
_BAM_EXT = (".bam", ".cram")


class UnknownTrackFormatError(ValueError):
    pass


def detect_format(path: str) -> tuple[TrackKind, str]:
    """Return (kind, format) for a file path."""
    p = os.fspath(path).lower()
    base = p[:-3] if p.endswith(".gz") else p
    if base.endswith(_CONTINUOUS_EXT):
        return TrackKind.CONTINUOUS, "bigwig"
    if base.endswith(_BEDGRAPH_EXT):
        return TrackKind.CONTINUOUS, "bedgraph"
    if base.endswith(_BAM_EXT):
        return TrackKind.CONTINUOUS, "bam"
    if base.endswith(_INTERVAL_EXT):
        return TrackKind.INTERVAL, detect_bed_format(base)
    raise UnknownTrackFormatError(
        f"{path}: unrecognized track format. Supported: bigWig (.bw/.bigWig), "
        "bedGraph (.bedGraph/.bg/.bdg), BED (.bed), narrowPeak, broadPeak, BAM/CRAM; "
        ".gz accepted for text formats."
    )


def open_track(
    path: str,
    *,
    name: str | None = None,
    kind: str | None = None,
    format: str | None = None,
    **kwargs,
) -> Track:
    """Open any supported track.

    ``kind``/``format`` override sniffing — e.g. ``kind='point'`` loads a
    BED file of modification sites as a PointTrack instead of intervals.
    """
    if format is None or kind is None:
        det_kind, det_format = detect_format(path)
        kind = kind or det_kind.value
        format = format or det_format
    kind = TrackKind(kind)
    if kind is TrackKind.POINT:
        return open_points(path, name=name, **kwargs)
    if format == "bigwig":
        return open_bigwig(path, name=name, **kwargs)
    if format == "bedgraph":
        return open_bedgraph(path, name=name, **kwargs)
    if format == "bam":
        from .bam import open_bam_coverage

        return open_bam_coverage(path, name=name, **kwargs)
    if format in ("bed", "narrowpeak", "broadpeak"):
        return open_intervals(path, name=name, format=format, **kwargs)
    raise UnknownTrackFormatError(f"{path}: unsupported format {format!r}")
