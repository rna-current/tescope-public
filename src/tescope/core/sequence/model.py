"""Consensus-sequence model: sequence + provenance + resolution status.

Scientific safety contract (see docs/sequence-browser.md):

- A sequence is displayed as coordinate-exact ONLY when its compatibility
  with the analysis' RepeatMasker-derived consensus model has been
  positively established (exact length match — the same non-inferential
  gate the feature-annotation layer uses).
- TEscope never invents bases, fills gaps, concatenates database models,
  or remaps an incompatible sequence onto the analysis axis.
- The sequence layer is display/metadata only: it takes no part in
  projection, aggregation, filtering or normalization.

Coordinates: sequences are stored 0-based (index i = consensus position
[i, i+1)); user-facing base positions are 1-based inclusive (see
:mod:`tescope.core.sequence.coords`).
"""

from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass, field

# resolution statuses (superset of the feature-layer Compatibility enum:
# sequence sources can also be entirely absent)
STATUS_VERIFIED = "verified"  # positively matched to the analysis model
STATUS_UNVERIFIED = "unverified"  # compatibility could not be established
STATUS_INCOMPATIBLE = "incompatible"  # positively determined NOT to match
STATUS_UNAVAILABLE = "unavailable"  # no source / no family / offline+uncached

SEQUENCE_STATUSES = (
    STATUS_VERIFIED,
    STATUS_UNVERIFIED,
    STATUS_INCOMPATIBLE,
    STATUS_UNAVAILABLE,
)

# IUPAC nucleotide codes (uppercase). Repeat consensus libraries commonly
# contain ambiguity codes; anything outside this set is rejected.
IUPAC_NUCLEOTIDES = frozenset("ACGTUNRYSWKMBDHV")

_COMPLEMENT = str.maketrans(
    "ACGTUNRYSWKMBDHV",
    "TGCAANYRSWMKVHDB",
)


class SequenceError(ValueError):
    """User-facing sequence-layer failure with an actionable message."""


def clean_sequence(raw: str, where: str = "sequence") -> str:
    """Uppercase and validate a nucleotide sequence (IUPAC alphabet).

    Whitespace is removed; any other non-IUPAC character fails with its
    1-based position — TEscope never silently drops or substitutes bases.
    """
    seq = "".join(raw.split()).upper()
    for i, ch in enumerate(seq):
        if ch not in IUPAC_NUCLEOTIDES:
            raise SequenceError(
                f"{where}: invalid character {ch!r} at base {i + 1} "
                "(expected IUPAC nucleotide codes)"
            )
    return seq


def sequence_checksum(seq: str) -> str:
    """sha256 over the uppercase sequence (identity for cache/provenance)."""
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class SequenceProvenance:
    """Where the consensus sequence came from — retained end to end."""

    provider: str = ""  # "Dfam" | "user-provided library" | "user-provided FASTA"
    family: str = ""  # the analysis' RepeatMasker repeat name
    record_name: str = ""  # source's own record/header name
    accession: str = ""  # e.g. "DF000000225" (Dfam)
    source_version: str = ""  # e.g. "Dfam 4.0" or user-entered library name
    source_release_date: str = ""
    retrieved_at: str = ""  # ISO timestamp (remote sources)
    api_url: str = ""  # endpoint (remote sources)
    source_path: str = ""  # local path (user-provided sources)
    sequence_length: int | None = None
    checksum: str = ""  # sha256 of the uppercase sequence
    expected_model_length: int | None = None  # TEscope consensus model (bp)
    matched_by: str = ""  # "name" | "alias:Repbase" | "header" | "only record"
    compatibility: str = STATUS_UNVERIFIED
    compatibility_detail: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> SequenceProvenance:
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        return cls(**known)


@dataclass
class ConsensusSequence:
    """A validated consensus sequence with its provenance."""

    sequence: str  # uppercase IUPAC, 0-based indexing
    provenance: SequenceProvenance = field(default_factory=SequenceProvenance)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def to_dict(self) -> dict:
        return {"sequence": self.sequence, "provenance": self.provenance.to_dict()}

    @classmethod
    def from_dict(cls, d: dict) -> ConsensusSequence:
        return cls(
            sequence=str(d.get("sequence", "")),
            provenance=SequenceProvenance.from_dict(d.get("provenance") or {}),
        )


@dataclass
class SequenceResolution:
    """Outcome of resolving a consensus sequence for one analysis.

    ``sequence`` is attached ONLY when the status is VERIFIED — an
    incompatible or unverified candidate is never handed to display code,
    so coordinate-exact functions cannot show wrong bases by accident.
    Provenance is kept in every case (it records why).
    """

    status: str = STATUS_UNAVAILABLE
    detail: str = ""
    sequence: ConsensusSequence | None = None
    provenance: SequenceProvenance | None = None

    @property
    def verified(self) -> bool:
        return self.status == STATUS_VERIFIED and self.sequence is not None

    def to_dict(self) -> dict:
        return {
            "status": self.status,
            "detail": self.detail,
            "sequence": self.sequence.to_dict() if self.sequence else None,
            "provenance": self.provenance.to_dict() if self.provenance else None,
        }

    @classmethod
    def from_dict(cls, d: dict) -> SequenceResolution:
        seq = ConsensusSequence.from_dict(d["sequence"]) if d.get("sequence") else None
        prov = (
            SequenceProvenance.from_dict(d["provenance"]) if d.get("provenance") else None
        )
        return cls(
            status=str(d.get("status", STATUS_UNAVAILABLE)),
            detail=str(d.get("detail", "")),
            sequence=seq,
            provenance=prov,
        )
