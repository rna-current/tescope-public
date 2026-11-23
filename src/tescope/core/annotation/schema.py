"""Canonical TEscope annotation schema.

One row per RepeatMasker hit. All coordinates are stored 0-based,
half-open ([start, end)), regardless of the source format:

- ``genome_start`` / ``genome_end``: position of the hit on the genome.
- ``cons_start`` / ``cons_end``: position of the hit on the repeat
  consensus/model, always expressed in the consensus's own 5'->3'
  orientation (RepeatMasker already reports consensus coordinates this
  way for both strands).
- ``cons_left``: bases of consensus remaining 3' of ``cons_end``
  (RepeatMasker's parenthesized "(left)" field). ``-1`` when unknown.

The consensus length of a hit is therefore ``cons_end + cons_left``.
Different hits of the same repeat name can disagree slightly (RepeatMasker
version drift); per-repeat consensus length is taken as the maximum over
hits (see :func:`consensus_lengths`).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone

import numpy as np
import pandas as pd

# Column name -> dtype of the canonical annotation frame.
CANONICAL_COLUMNS: dict[str, str] = {
    "chrom": "category",
    "genome_start": "int64",
    "genome_end": "int64",
    "strand": "category",  # '+' or '-'
    "repeat_name": "category",
    "repeat_class": "category",
    "repeat_family": "category",
    "sw_score": "float64",  # NaN if absent
    "pct_div": "float64",
    "pct_del": "float64",
    "pct_ins": "float64",
    "cons_start": "int64",
    "cons_end": "int64",
    "cons_left": "int64",  # -1 if unknown
    "rm_id": "int64",  # RepeatMasker linkage ID, -1 if absent
    "is_secondary": "bool",  # trailing '*' in .out (overlaps a higher-scoring hit)
}


@dataclass
class AnnotationProvenance:
    """Where an annotation set came from. Stored alongside every cache."""

    assembly: str
    source: str  # e.g. "ucsc-goldenpath", "user-file"
    source_url: str = ""
    source_path: str = ""
    source_format: str = ""  # "repeatmasker-out" | "ucsc-rmsk"
    retrieved_at: str = ""  # ISO timestamp
    source_version: str = ""  # RepeatMasker/library version if known
    md5: str = ""
    notes: str = ""
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> AnnotationProvenance:
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        return cls(**known)

    @staticmethod
    def now() -> str:
        return datetime.now(timezone.utc).isoformat(timespec="seconds")


class AnnotationSchemaError(ValueError):
    """Raised when a frame does not conform to the canonical schema."""


def empty_frame() -> pd.DataFrame:
    """Return an empty canonical annotation frame."""
    data = {name: pd.Series(dtype=dt) for name, dt in CANONICAL_COLUMNS.items()}
    return pd.DataFrame(data)


def finalize_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Cast to canonical dtypes, order columns, and validate invariants."""
    missing = [c for c in CANONICAL_COLUMNS if c not in df.columns]
    if missing:
        raise AnnotationSchemaError(f"annotation frame missing columns: {missing}")
    out = df[list(CANONICAL_COLUMNS)].copy()
    for name, dt in CANONICAL_COLUMNS.items():
        out[name] = out[name].astype(dt)
    validate_frame(out)
    return out.reset_index(drop=True)


def validate_frame(df: pd.DataFrame) -> None:
    """Cheap invariant checks; raises AnnotationSchemaError on violation."""
    if len(df) == 0:
        return
    if not set(df["strand"].unique()).issubset({"+", "-"}):
        bad = sorted(set(df["strand"].unique()) - {"+", "-"})
        raise AnnotationSchemaError(f"invalid strand values: {bad}")
    if (df["genome_end"] <= df["genome_start"]).any():
        n = int((df["genome_end"] <= df["genome_start"]).sum())
        raise AnnotationSchemaError(f"{n} hits with genome_end <= genome_start")
    if (df["cons_end"] < df["cons_start"]).any():
        n = int((df["cons_end"] < df["cons_start"]).sum())
        raise AnnotationSchemaError(f"{n} hits with cons_end < cons_start")
    if (df["cons_start"] < 0).any():
        raise AnnotationSchemaError("negative cons_start")


def consensus_lengths(df: pd.DataFrame) -> pd.Series:
    """Per-repeat-name consensus/model length.

    Length of a hit's model is ``cons_end + cons_left`` (unknown ``cons_left``
    contributes only ``cons_end``). Real annotations carry rare artifact hits
    with inflated consensus coordinates (e.g. L1HS hits implying a 7.8 kb
    model when the true model is ~6.1 kb), so the **mode** (most frequent
    value; largest on ties) across hits is used, not the maximum. Hits whose
    coordinates exceed the model length are excluded at analysis time with a
    warning.
    """
    hit_len = df["cons_end"] + df["cons_left"].clip(lower=0)

    def _mode_max(s: pd.Series) -> int:
        m = s.mode()
        return int(m.max()) if len(m) else 0

    return hit_len.groupby(df["repeat_name"], observed=True).agg(_mode_max)


def repeat_hierarchy(df: pd.DataFrame) -> pd.DataFrame:
    """Class -> family -> name hierarchy with copy counts."""
    g = (
        df.groupby(["repeat_class", "repeat_family", "repeat_name"], observed=True)
        .size()
        .rename("n_copies")
        .reset_index()
    )
    return g.sort_values(["repeat_class", "repeat_family", "repeat_name"]).reset_index(drop=True)


def add_derived_columns(df: pd.DataFrame, cons_len: pd.Series | None = None) -> pd.DataFrame:
    """Attach derived per-hit columns used by filtering and projection.

    Adds:
      - ``genome_length``
      - ``cons_span``: consensus bases covered by the hit
      - ``cons_model_length``: model length for the hit's repeat name
      - ``consensus_coverage``: cons_span / cons_model_length
    """
    out = df.copy()
    out["genome_length"] = out["genome_end"] - out["genome_start"]
    out["cons_span"] = out["cons_end"] - out["cons_start"]
    if cons_len is None:
        cons_len = consensus_lengths(df)
    model_len = out["repeat_name"].map(cons_len).astype("float64")
    out["cons_model_length"] = model_len
    with np.errstate(divide="ignore", invalid="ignore"):
        cov = out["cons_span"] / model_len
    out["consensus_coverage"] = cov.clip(upper=1.0)
    return out


def drop_inconsistent_hits(df: pd.DataFrame, cons_length: int) -> tuple[pd.DataFrame, int]:
    """Exclude hits whose consensus coordinates exceed the model length.

    Real annotations carry rare artifact hits beyond the mode-derived model
    length; projecting them would pile signal into the last bin. Shared by
    the analysis orchestrator and the GUI's live filter counts so both
    always report the same copy set.
    """
    bad = df["cons_end"] > cons_length
    n = int(bad.sum())
    return (df[~bad], n) if n else (df, 0)
