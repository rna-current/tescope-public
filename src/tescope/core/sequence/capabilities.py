"""Sequence capabilities: the one place answering "what sequence
functionality is genuinely available right now?".

Extends the centralized capability philosophy of
:mod:`tescope.core.tracks.capabilities` to the sequence dimension. GUI
controls and CLI checks consult this — conditions are never scattered
through widget code. Every unavailable capability carries a
human-readable reason used as the tooltip/error text.

Sequence capability is ORTHOGONAL to signal-track capability: a bigWig
analysis can have a verified sequence track, and a BAM analysis can have
none. Coordinate-exact functions (base display, copy, FASTA export,
search) activate only for a VERIFIED sequence.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .model import (
    STATUS_INCOMPATIBLE,
    STATUS_UNAVAILABLE,
    STATUS_VERIFIED,
    SequenceResolution,
)

_NOT_CONFIGURED = (
    "No consensus sequence source is configured. Choose Automatic (Dfam) "
    "or a user-provided FASTA/library in the Sequence panel."
)
_INCOMPATIBLE = (
    "The available sequence does not match this analysis' consensus "
    "coordinate model, so bases cannot be shown at exact coordinates. "
)


@dataclass(frozen=True)
class SequenceCapabilities:
    """What the current sequence resolution genuinely supports."""

    available: bool = False  # a verified sequence is attached
    verified: bool = False  # compatibility positively established
    base_display: bool = False  # base-level letters at exact coordinates
    copy_sequence: bool = False
    export_fasta: bool = False
    search: bool = False
    status: str = STATUS_UNAVAILABLE
    # capability name -> reason it is unavailable (tooltip text)
    reasons: dict = field(default_factory=dict)

    def why_not(self, capability: str) -> str:
        return self.reasons.get(
            capability, self.reasons.get("all", "Not available.")
        )


_EXACT_CAPS = ("base_display", "copy_sequence", "export_fasta", "search")


def sequence_capabilities(
    res: SequenceResolution | None,
) -> SequenceCapabilities:
    """Capabilities for one resolution (None = no resolution attempted)."""
    if res is None:
        reason = _NOT_CONFIGURED
        return SequenceCapabilities(
            status=STATUS_UNAVAILABLE,
            reasons={"all": reason, **{c: reason for c in _EXACT_CAPS}},
        )
    if res.status == STATUS_VERIFIED and res.sequence is not None:
        return SequenceCapabilities(
            available=True,
            verified=True,
            base_display=True,
            copy_sequence=True,
            export_fasta=True,
            search=True,
            status=res.status,
        )
    if res.status == STATUS_INCOMPATIBLE:
        reason = _INCOMPATIBLE + (res.detail or "")
    else:
        reason = res.detail or _NOT_CONFIGURED
    return SequenceCapabilities(
        status=res.status,
        reasons={"all": reason, **{c: reason for c in _EXACT_CAPS}},
    )
