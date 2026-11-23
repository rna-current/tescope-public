"""Upstream per-locus estimate import (tescope-locus-estimates v1).

Scientific contract (docs/feature-import.md, docs/dev/
SCIENTIFIC_CONTRACT_locus_features.md):

- An estimate table supports locus-level TABLES and comparisons only.
  It is never drawn as positional signal and never combined on an
  axis with coverage units.
- ``total`` may already include unique evidence
  (``total_includes_unique: true`` — e.g. Telescope ``final_count``
  with ``unique_count``); TEscope never adds the two. A missing
  ``unique`` column is *unavailable*, not zero.
- Native score semantics are recorded verbatim and never normalized
  across tools. Selecting a locus subset never re-fits or reassigns
  anything.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd

_HEADER_RE = re.compile(r"^#\s*([A-Za-z_]+)\s*:\s*(.*?)\s*$")
_MAGIC = "tescope-locus-estimates"
_REQUIRED_DIRECTIVES = ("assay", "units", "total_includes_unique")
_LOCUS_ID_RE = re.compile(r"^[^:\s]+:\d+-\d+\([+-]\)$")


class EstimateError(ValueError):
    """User-facing estimate-import failure."""


@dataclass
class EstimateMatchReport:
    n_rows: int = 0
    n_matched: int = 0
    unmatched_ids: list[str] = field(default_factory=list)  # first 20
    n_unmatched: int = 0
    assembly_checked: bool = False

    def summary(self) -> str:
        s = f"{self.n_matched}/{self.n_rows} estimate rows matched analyzed loci"
        if self.n_unmatched:
            s += f"; {self.n_unmatched} unmatched (reported, not dropped silently)"
        return s


@dataclass
class LocusEstimates:
    path: str
    tool: str  # recorded verbatim; "" = unknown
    assay: str
    units: str
    total_includes_unique: bool | None  # None = unknown (declared)
    assembly: str  # "" = not declared
    table: pd.DataFrame = field(default_factory=pd.DataFrame)
    # columns: locus_id (index), total, unique (optional), score (optional)

    @property
    def has_unique(self) -> bool:
        return "unique" in self.table.columns

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "tool": self.tool,
            "assay": self.assay,
            "units": self.units,
            "total_includes_unique": self.total_includes_unique,
            "assembly": self.assembly,
        }


def load_estimates(path: str) -> LocusEstimates:
    directives: dict[str, str] = {}
    header_line = None
    data_lines: list[str] = []
    with open(path, encoding="utf-8") as fh:
        first = fh.readline()
        if _MAGIC not in first:
            raise EstimateError(
                "not a tescope-locus-estimates file: first line must contain "
                f"'{_MAGIC} v1' (see docs/feature-import.md)"
            )
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            m = _HEADER_RE.match(line)
            if m and header_line is None:
                directives[m.group(1).lower()] = m.group(2)
                continue
            if header_line is None:
                header_line = line
                continue
            data_lines.append(line)
    for key in _REQUIRED_DIRECTIVES:
        if key not in directives:
            raise EstimateError(
                f"missing required header directive '# {key}:' — unknown "
                "semantics are not guessed (see docs/feature-import.md)"
            )
    tiu_raw = directives["total_includes_unique"].strip().lower()
    if tiu_raw in ("true", "yes"):
        tiu: bool | None = True
    elif tiu_raw in ("false", "no"):
        tiu = False
    elif tiu_raw == "unknown":
        tiu = None
    else:
        raise EstimateError(
            "total_includes_unique must be true/false/unknown, "
            f"got {tiu_raw!r}"
        )
    if directives["units"].strip().lower() in ("aa", "amino-acids", "aminoacids"):
        raise EstimateError("amino-acid units are not supported")
    if header_line is None:
        raise EstimateError("no column header line found")
    columns = header_line.split("\t")
    if "locus_id" not in columns or "total" not in columns:
        raise EstimateError(
            "column header must contain at least 'locus_id' and 'total', "
            f"got {columns}"
        )
    rows = [line.split("\t") for line in data_lines]
    bad = [i for i, r in enumerate(rows) if len(r) != len(columns)]
    if bad:
        raise EstimateError(
            f"row {bad[0] + 1}: expected {len(columns)} tab-separated "
            "fields (malformed row)"
        )
    df = pd.DataFrame(rows, columns=columns)
    bad_ids = [x for x in df["locus_id"] if not _LOCUS_ID_RE.match(str(x))]
    if bad_ids:
        raise EstimateError(
            f"invalid locus_id {bad_ids[0]!r}: expected "
            "'chrom:start-end(strand)' (0-based half-open; export the "
            "analyzed IDs from the Loci panel)"
        )
    dup = df["locus_id"][df["locus_id"].duplicated()]
    if len(dup):
        raise EstimateError(f"duplicate locus_id {dup.iloc[0]!r}")
    for col in ("total", "unique", "score"):
        if col in df.columns:
            try:
                df[col] = pd.to_numeric(df[col])
            except (ValueError, TypeError) as e:
                raise EstimateError(f"column '{col}' is not numeric: {e}") from None
    if "unique" in df.columns and tiu is True:
        over = df[df["unique"] > df["total"]]
        if len(over):
            raise EstimateError(
                f"locus {over['locus_id'].iloc[0]!r}: unique > total although "
                "total_includes_unique is true — inconsistent decomposition"
            )
    df = df.set_index("locus_id")
    return LocusEstimates(
        path=path,
        tool=directives.get("tool", ""),
        assay=directives["assay"],
        units=directives["units"],
        total_includes_unique=tiu,
        assembly=directives.get("assembly", ""),
        table=df,
    )


def match_estimates(
    est: LocusEstimates, locus_table: pd.DataFrame, assembly: str
) -> tuple[pd.DataFrame, EstimateMatchReport]:
    """Join estimate columns onto the locus table as est_* columns.

    Returns (augmented table, report). Unmatched estimate rows are
    reported; analyzed loci without an estimate get NaN (unavailable,
    not zero).
    """
    if est.assembly and est.assembly != assembly:
        raise EstimateError(
            f"estimate table declares assembly {est.assembly!r} but the "
            f"analysis uses {assembly!r}; explicit reconversion is required "
            "(no silent aliasing)"
        )
    report = EstimateMatchReport(
        n_rows=len(est.table), assembly_checked=bool(est.assembly)
    )
    known = set(locus_table["locus_id"])
    unmatched = [i for i in est.table.index if i not in known]
    report.n_unmatched = len(unmatched)
    report.unmatched_ids = [str(x) for x in unmatched[:20]]
    report.n_matched = report.n_rows - report.n_unmatched

    out = locus_table.copy()
    joined = est.table.reindex(out["locus_id"])
    out["est_total"] = joined["total"].to_numpy()
    if "unique" in joined.columns:
        out["est_unique"] = joined["unique"].to_numpy()
    if "score" in joined.columns:
        out["est_score"] = joined["score"].to_numpy()
    return out, report
