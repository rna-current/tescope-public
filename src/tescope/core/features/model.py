"""Consensus-feature annotation model.

Scientific safety contract (see docs/annotations.md):

- TEscope NEVER invents, infers, or hard-codes biological feature names or
  coordinates. Features come verbatim from an authoritative source record
  (with provenance) or are explicitly user-defined.
- Source terminology is preserved: a Dfam ``protein_match`` stays a
  protein match — it is never promoted to "ORF"/"gene".
- Features are a metadata/visualization layer. They are attached to
  results but take no part in projection, aggregation, or normalization.

Coordinates: internally 0-based half-open on the TEscope consensus axis;
the original source coordinates and convention are preserved in
provenance.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum


class Compatibility(str, Enum):
    """Whether external feature coordinates fit the analysis consensus.

    - VERIFIED: the source consensus was positively matched to the
      annotation-derived model (exact length match).
    - INCOMPATIBLE: the source consensus was positively determined NOT to
      match (e.g. different length). Features must not be displayed on
      the analysis axis.
    - UNVERIFIED: compatibility could not be established (missing data).
      Treated like incompatible for display unless the user opts in
      explicitly; TEscope itself never displays unverified coordinates.
    """

    VERIFIED = "verified"
    INCOMPATIBLE = "incompatible"
    UNVERIFIED = "unverified"


# origin of a feature record — never conflated
SOURCE_DATABASE = "database"
SOURCE_USER = "user-defined"


@dataclass
class ConsensusFeature:
    """One feature interval on the consensus axis."""

    name: str  # display label, verbatim from source or user
    start: int  # 0-based half-open on the TEscope consensus
    end: int
    origin: str  # SOURCE_DATABASE | SOURCE_USER
    type: str = ""  # source-provided type (e.g. "protein_match", "cds"); "" for user
    note: str = ""  # free text (user) or source description
    # source-record details (database features only)
    source_attributes: dict = field(default_factory=dict)

    def validate(self, cons_length: int | None = None) -> None:
        if not self.name:
            raise ValueError("feature name must not be empty")
        if self.end <= self.start:
            raise ValueError(
                f"feature '{self.name}': end ({self.end}) must be > start ({self.start})"
            )
        if self.start < 0:
            raise ValueError(f"feature '{self.name}': negative start")
        if cons_length is not None and self.end > cons_length:
            raise ValueError(
                f"feature '{self.name}': end ({self.end}) exceeds consensus "
                f"length ({cons_length})"
            )

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> ConsensusFeature:
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        return cls(**known)


@dataclass
class FeatureProvenance:
    """Where database features came from — retained end to end."""

    source: str = ""  # e.g. "Dfam"
    accession: str = ""  # e.g. "DF000004146"
    record_name: str = ""  # source's own family name
    record_version: str = ""  # per-record version where provided
    database_version: str = ""  # e.g. "Dfam 4.0"
    database_release_date: str = ""
    retrieved_at: str = ""  # ISO timestamp (UTC)
    api_url: str = ""
    source_coordinate_convention: str = ""  # e.g. "1-based inclusive (model_start/end_pos)"
    source_consensus_length: int | None = None
    matched_by: str = ""  # "name" | "alias:Repbase" | ...
    compatibility: str = Compatibility.UNVERIFIED.value
    compatibility_detail: str = ""  # why accepted/rejected

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> FeatureProvenance:
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        return cls(**known)


@dataclass
class FeatureSet:
    """Features attached to one analysis (one repeat, one consensus axis)."""

    repeat_name: str
    cons_length: int
    database_features: list[ConsensusFeature] = field(default_factory=list)
    user_features: list[ConsensusFeature] = field(default_factory=list)
    provenance: FeatureProvenance | None = None
    status: str = ""  # human-readable annotation state for QC/GUI

    @property
    def n_total(self) -> int:
        return len(self.database_features) + len(self.user_features)

    def all_features(self) -> list[ConsensusFeature]:
        return list(self.database_features) + list(self.user_features)

    def validate(self) -> None:
        for f in self.all_features():
            f.validate(self.cons_length)
        for f in self.database_features:
            if f.origin != SOURCE_DATABASE:
                raise ValueError(f"feature '{f.name}' in database lane is not database-origin")
        for f in self.user_features:
            if f.origin != SOURCE_USER:
                raise ValueError(f"feature '{f.name}' in user lane is not user-origin")

    def to_dict(self) -> dict:
        return {
            "repeat_name": self.repeat_name,
            "cons_length": self.cons_length,
            "database_features": [f.to_dict() for f in self.database_features],
            "user_features": [f.to_dict() for f in self.user_features],
            "provenance": self.provenance.to_dict() if self.provenance else None,
            "status": self.status,
        }

    @classmethod
    def from_dict(cls, d: dict) -> FeatureSet:
        return cls(
            repeat_name=d["repeat_name"],
            cons_length=int(d["cons_length"]),
            database_features=[ConsensusFeature.from_dict(x) for x in d.get("database_features", [])],
            user_features=[ConsensusFeature.from_dict(x) for x in d.get("user_features", [])],
            provenance=FeatureProvenance.from_dict(d["provenance"]) if d.get("provenance") else None,
            status=d.get("status", ""),
        )


def verify_length_compatibility(
    expected_length: int, source_length: int | None
) -> tuple[Compatibility, str]:
    """Compatibility by exact consensus-length match.

    The TEscope model length is derived from the RepeatMasker annotation
    (mode of end+left); a source consensus of a different length cannot be
    assumed to share coordinates, so anything but an exact match is
    INCOMPATIBLE. Missing information is UNVERIFIED, never silently
    accepted.
    """
    if source_length is None or expected_length <= 0:
        return (
            Compatibility.UNVERIFIED,
            "source or annotation consensus length unavailable",
        )
    if int(source_length) == int(expected_length):
        return (
            Compatibility.VERIFIED,
            f"consensus length match ({expected_length} bp)",
        )
    return (
        Compatibility.INCOMPATIBLE,
        f"consensus length mismatch: annotation model {expected_length} bp "
        f"vs source consensus {source_length} bp",
    )
