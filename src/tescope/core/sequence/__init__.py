"""Consensus-sequence support.

Provides the sequence displayed by the consensus browser. The displayed
sequence must correspond to the SAME consensus coordinate model TEscope
projects onto — family-name equality never implies model equality, so
every candidate sequence passes the same exact-length compatibility gate
the feature-annotation subsystem uses. TEscope never invents, joins or
remaps sequences to force a fit.

Sources: Dfam (automatic, official API, cached), user-provided FASTA
libraries/files (read locally, never redistributed), or none.
"""

from .model import (
    STATUS_INCOMPATIBLE,
    STATUS_UNAVAILABLE,
    STATUS_UNVERIFIED,
    STATUS_VERIFIED,
    ConsensusSequence,
    SequenceProvenance,
    SequenceResolution,
)
from .resolve import resolve_sequence

__all__ = [
    "STATUS_INCOMPATIBLE",
    "STATUS_UNAVAILABLE",
    "STATUS_UNVERIFIED",
    "STATUS_VERIFIED",
    "ConsensusSequence",
    "SequenceProvenance",
    "SequenceResolution",
    "resolve_sequence",
]
