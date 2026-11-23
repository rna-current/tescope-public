"""TE copy filtering: genomic length, consensus coverage, completeness.

Filters operate on the canonical annotation frame (plus derived columns from
:func:`tescope.core.annotation.schema.add_derived_columns`) and are pure:
they return boolean masks / filtered frames and a step-by-step count report
so both GUI and CLI can show "N copies remain after each filter".

Completeness classification is geometry-only (where the hit sits on the
consensus model) and deliberately biology-agnostic; tolerances are
parameters, not hard-coded per-family rules.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

COMPLETENESS_CLASSES = (
    "full_length",
    "five_prime_truncated",
    "three_prime_truncated",
    "internal_fragment",
)


@dataclass
class FilterConfig:
    """User-facing filter settings (all optional; None/empty = no filter)."""

    min_genome_length: int | None = None
    min_consensus_coverage: float | None = None  # 0..1
    completeness: list[str] = field(default_factory=list)  # subset of COMPLETENESS_CLASSES
    # Tolerance for calling an end "reached": a hit is considered to reach the
    # consensus 5'/3' end if it comes within `end_tolerance_frac * model_length`
    # or `end_tolerance_bp` (whichever is larger) of that end.
    end_tolerance_frac: float = 0.05
    end_tolerance_bp: int = 20
    include_secondary: bool = True  # keep hits flagged '*' by RepeatMasker
    strands: list[str] = field(default_factory=list)  # e.g. ["+"]; empty = both

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> FilterConfig:
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        bad = set(known.get("completeness") or []) - set(COMPLETENESS_CLASSES)
        if bad:
            raise ValueError(
                f"unknown completeness classes {sorted(bad)}; valid: {COMPLETENESS_CLASSES}"
            )
        return cls(**known)


def classify_completeness(
    df: pd.DataFrame,
    *,
    end_tolerance_frac: float = 0.05,
    end_tolerance_bp: int = 20,
) -> pd.Series:
    """Classify each hit by which consensus ends it reaches.

    Requires derived columns (``cons_model_length``). Hits whose model
    length is unknown (NaN/0) are classified by consensus span only and
    fall back to ``internal_fragment``.
    """
    model_len = df["cons_model_length"].to_numpy(dtype="float64")
    tol = np.maximum(model_len * end_tolerance_frac, end_tolerance_bp)
    reaches_5p = df["cons_start"].to_numpy() <= tol
    with np.errstate(invalid="ignore"):
        reaches_3p = (model_len - df["cons_end"].to_numpy()) <= tol
    ok_model = np.isfinite(model_len) & (model_len > 0)
    reaches_5p = reaches_5p & ok_model
    reaches_3p = reaches_3p & ok_model

    out = np.full(len(df), "internal_fragment", dtype=object)
    out[reaches_5p & reaches_3p] = "full_length"
    out[reaches_5p & ~reaches_3p] = "three_prime_truncated"
    out[~reaches_5p & reaches_3p] = "five_prime_truncated"
    return pd.Series(pd.Categorical(out, categories=list(COMPLETENESS_CLASSES)), index=df.index)


@dataclass
class FilterStep:
    label: str
    n_remaining: int


@dataclass
class FilterReport:
    n_input: int
    steps: list[FilterStep] = field(default_factory=list)

    @property
    def n_remaining(self) -> int:
        return self.steps[-1].n_remaining if self.steps else self.n_input

    def to_rows(self) -> list[tuple[str, int]]:
        rows = [("All copies", self.n_input)]
        rows += [(s.label, s.n_remaining) for s in self.steps]
        return rows


def apply_filters(df: pd.DataFrame, cfg: FilterConfig) -> tuple[pd.DataFrame, FilterReport]:
    """Apply filters sequentially; returns (filtered frame, count report).

    ``df`` must carry the derived columns (``genome_length``,
    ``consensus_coverage``, ``cons_model_length``).
    """
    report = FilterReport(n_input=len(df))
    mask = pd.Series(True, index=df.index)

    if not cfg.include_secondary and "is_secondary" in df.columns:
        mask &= ~df["is_secondary"]
        report.steps.append(FilterStep("Primary hits only", int(mask.sum())))

    if cfg.strands:
        mask &= df["strand"].isin(cfg.strands)
        report.steps.append(FilterStep(f"Strand in {cfg.strands}", int(mask.sum())))

    if cfg.min_genome_length is not None:
        mask &= df["genome_length"] >= cfg.min_genome_length
        report.steps.append(
            FilterStep(f"Length >= {_fmt_bp(cfg.min_genome_length)}", int(mask.sum()))
        )

    if cfg.min_consensus_coverage is not None:
        cov_ok = df["consensus_coverage"] >= cfg.min_consensus_coverage
        cov_ok &= df["consensus_coverage"].notna()
        mask &= cov_ok
        pct = int(round(cfg.min_consensus_coverage * 100))
        report.steps.append(FilterStep(f"Consensus coverage >= {pct}%", int(mask.sum())))

    if cfg.completeness:
        cls = classify_completeness(
            df,
            end_tolerance_frac=cfg.end_tolerance_frac,
            end_tolerance_bp=cfg.end_tolerance_bp,
        )
        mask &= cls.isin(cfg.completeness)
        label = "Completeness: " + ", ".join(c.replace("_", " ") for c in cfg.completeness)
        report.steps.append(FilterStep(label, int(mask.sum())))

    return df[mask], report


def _fmt_bp(n: int) -> str:
    if n >= 1000 and n % 1000 == 0:
        return f"{n // 1000} kb"
    return f"{n} bp"
