"""Input-aware track capabilities.

One central place answering "what can this input actually support?", so
GUI controls and config validation never scatter file-extension checks.
Capabilities are probed from the ACTUAL file where practical (e.g. a BAM
is inspected for NH tags rather than assumed to have them); a bigWig can
never acquire mapping information, and TEscope refuses to pretend
otherwise.

Every unavailable capability carries a human-readable reason used as the
GUI tooltip and in config-validation errors.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

from .loader import UnknownTrackFormatError, detect_format

_PREAGG_REASON = (
    "Unavailable for pre-aggregated signal tracks. Mapping classification "
    "requires read-level alignment information such as BAM/CRAM."
)
_NO_NH_REASON = "Unavailable: this alignment file does not contain NH tags."
_NOT_CONTINUOUS_SMOOTH = (
    "Smoothing applies to continuous signal profiles; this track produces "
    "interval/point occurrence output."
)
_NO_PYSAM = "BAM/CRAM support requires the pysam extra (unavailable on Windows)."

# how many reads to inspect when probing BAM tags
_PROBE_READS = 200


@dataclass(frozen=True)
class TrackCapabilities:
    """What one input file genuinely supports."""

    kind: str  # continuous | interval | point
    format: str
    continuous_signal: bool = False
    interval_signal: bool = False
    point_signal: bool = False
    per_copy_projection: bool = False  # continuous tracks get a CopyBinMatrix
    smoothing: bool = False
    copy_distribution: bool = False  # per-copy quantile bands
    heatmap: bool = False
    occurrence: bool = False
    read_level: bool = False
    mapq: bool = False
    nh_tag: bool = False
    mapping_classification: bool = False  # any policy usable
    # capability name -> reason it is unavailable (tooltip text)
    reasons: dict = field(default_factory=dict)

    def why_not(self, capability: str) -> str:
        return self.reasons.get(capability, "Not supported by this input type.")


def _bam_probe(path: str) -> tuple[bool, bool]:
    """(openable, has_NH) by inspecting the first reads of the actual file."""
    try:
        import pysam
    except ImportError:
        return False, False
    mode = "rc" if str(path).lower().endswith(".cram") else "rb"
    try:
        with pysam.AlignmentFile(path, mode, check_sq=False) as af:
            for i, read in enumerate(af.fetch(until_eof=True)):
                if read.has_tag("NH"):
                    return True, True
                if i >= _PROBE_READS:
                    break
        return True, False
    except (OSError, ValueError):
        return False, False


def probe_capabilities(path: str, kind: str | None = None,
                       fmt: str | None = None) -> TrackCapabilities:
    """Capabilities of one input file (probing content where practical)."""
    if kind is None or fmt is None:
        det_kind, det_fmt = detect_format(path)
        kind = kind or det_kind.value
        fmt = fmt or det_fmt

    reasons: dict = {}
    if kind == "continuous" and fmt == "bam":
        openable, has_nh = (
            _bam_probe(path) if os.path.exists(path) else (True, False)
        )
        try:
            import pysam  # noqa: F401
            pysam_ok = True
        except ImportError:
            pysam_ok = False
        if not pysam_ok:
            reasons["mapping_classification"] = _NO_PYSAM
            reasons["mapq"] = _NO_PYSAM
            reasons["nh_tag"] = _NO_PYSAM
        elif not has_nh:
            reasons["nh_tag"] = _NO_NH_REASON
        return TrackCapabilities(
            kind=kind, format=fmt,
            continuous_signal=True, per_copy_projection=True, smoothing=True,
            copy_distribution=True, heatmap=True,
            read_level=pysam_ok, mapq=pysam_ok,
            nh_tag=pysam_ok and has_nh,
            mapping_classification=pysam_ok,
            reasons=reasons,
        )
    if kind == "continuous":
        reasons["mapping_classification"] = _PREAGG_REASON
        reasons["mapq"] = _PREAGG_REASON
        reasons["nh_tag"] = _PREAGG_REASON
        reasons["read_level"] = _PREAGG_REASON
        return TrackCapabilities(
            kind=kind, format=fmt,
            continuous_signal=True, per_copy_projection=True, smoothing=True,
            copy_distribution=True, heatmap=True, reasons=reasons,
        )
    if kind == "interval":
        reasons["smoothing"] = _NOT_CONTINUOUS_SMOOTH
        reasons["mapping_classification"] = _PREAGG_REASON
        return TrackCapabilities(
            kind=kind, format=fmt, interval_signal=True, occurrence=True,
            reasons=reasons,
        )
    # point
    reasons["smoothing"] = _NOT_CONTINUOUS_SMOOTH
    reasons["mapping_classification"] = _PREAGG_REASON
    return TrackCapabilities(
        kind=kind, format=fmt, point_signal=True, occurrence=True,
        reasons=reasons,
    )


def probe_config_track(tcfg) -> TrackCapabilities:
    """Capabilities for a TrackConfig (uses declared kind/format overrides)."""
    try:
        return probe_capabilities(
            tcfg.file, kind=tcfg.kind or None, fmt=tcfg.format or None
        )
    except UnknownTrackFormatError:
        return TrackCapabilities(kind="", format="", reasons={
            "all": "Unrecognized track format."})
