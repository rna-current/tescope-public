"""Reproducible locus-selection files (YAML).

A selection is a list of stable locus IDs plus the provenance needed
to interpret them (assembly, repeat, annotation identity). Importing
into a different context reports every non-matching ID instead of
silently adapting.
"""

from __future__ import annotations

import yaml

from ..analysis import AnalysisResult
from .table import locus_ids

SCHEMA = "tescope-locus-selection/1"


class SelectionError(ValueError):
    pass


def export_selection(
    result: AnalysisResult, selected: list[str], path: str
) -> None:
    known = set(locus_ids(result.copies))
    unknown = [s for s in selected if s not in known]
    if unknown:
        raise SelectionError(
            f"refusing to export: {len(unknown)} ID(s) not in the current "
            f"eligible population (first: {unknown[0]!r})"
        )
    doc = {
        "schema": SCHEMA,
        "assembly": result.config.assembly,
        "repeat": result.repeat_name,
        "coordinates": "0-based half-open",
        "annotation_provenance": dict(result.provenance or {}),
        "n_eligible_at_export": int(len(result.copies)),
        "locus_ids": list(selected),
    }
    with open(path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, allow_unicode=True)


def import_selection(
    path: str, result: AnalysisResult
) -> tuple[list[str], list[str], str]:
    """Return (valid_ids, invalid_ids, context_warning)."""
    with open(path, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    if not isinstance(doc, dict) or doc.get("schema") != SCHEMA:
        raise SelectionError(
            f"not a {SCHEMA} file (see the Loci panel's export action)"
        )
    ids = doc.get("locus_ids")
    if not isinstance(ids, list) or not all(isinstance(x, str) for x in ids):
        raise SelectionError("locus_ids must be a list of strings")
    warning = ""
    if doc.get("assembly") and doc["assembly"] != result.config.assembly:
        raise SelectionError(
            f"selection was exported for assembly {doc['assembly']!r}; "
            f"current analysis uses {result.config.assembly!r}"
        )
    if doc.get("repeat") and doc["repeat"] != result.repeat_name:
        warning = (
            f"selection was exported for repeat {doc['repeat']!r} "
            f"(current: {result.repeat_name!r})"
        )
    known = set(locus_ids(result.copies))
    valid = [s for s in ids if s in known]
    invalid = [s for s in ids if s not in known]
    return valid, invalid, warning
