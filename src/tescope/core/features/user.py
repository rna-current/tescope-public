"""User-defined consensus features.

Users may annotate the consensus axis themselves. Such features are
always labeled ``user-defined`` and are never presented as database
annotations.

File formats (coordinates are **0-based half-open** on the TEscope
consensus axis — the same convention as every TEscope coordinate; this is
documented in docs/annotations.md and validated against the consensus
length at attach time):

TSV (header required; ``label``/``note`` optional)::

    name\tstart\tend\tlabel\tnote
    my_region\t100\t900\tregulatory\tfrom our 2024 mapping

YAML/JSON (a list of mappings)::

    - name: my_region
      start: 100
      end: 900
      label: regulatory     # optional
      note: from our 2024 mapping   # optional
"""

from __future__ import annotations

import os
from pathlib import Path

import yaml

from .model import SOURCE_USER, ConsensusFeature

_REQUIRED = ("name", "start", "end")
_OPTIONAL = ("label", "note")


class UserFeatureError(ValueError):
    pass


def feature_from_record(rec: dict, where: str = "") -> ConsensusFeature:
    missing = [k for k in _REQUIRED if k not in rec or rec[k] in (None, "")]
    if missing:
        raise UserFeatureError(f"{where}: missing required field(s) {missing}")
    unknown = set(rec) - set(_REQUIRED) - set(_OPTIONAL)
    if unknown:
        raise UserFeatureError(f"{where}: unknown field(s) {sorted(unknown)}")
    try:
        start = int(rec["start"])
        end = int(rec["end"])
    except (TypeError, ValueError):
        raise UserFeatureError(
            f"{where}: start/end must be integers (0-based half-open)"
        ) from None
    f = ConsensusFeature(
        name=str(rec["name"]),
        start=start,
        end=end,
        origin=SOURCE_USER,
        type=str(rec.get("label", "") or ""),
        note=str(rec.get("note", "") or ""),
    )
    try:
        f.validate()  # bounds vs consensus length checked at attach time
    except ValueError as e:
        raise UserFeatureError(f"{where}: {e}") from e
    return f


def parse_user_features_file(path: str | os.PathLike) -> list[ConsensusFeature]:
    p = Path(path)
    if not p.exists():
        raise UserFeatureError(f"user feature file not found: {p}")
    if p.suffix.lower() in (".yaml", ".yml", ".json"):
        try:
            data = yaml.safe_load(p.read_text(encoding="utf-8"))
        except yaml.YAMLError as e:
            raise UserFeatureError(f"{p}: cannot parse: {e}") from e
        if not isinstance(data, list):
            raise UserFeatureError(f"{p}: expected a list of feature mappings")
        return [
            feature_from_record(rec, where=f"{p.name}[{i}]")
            for i, rec in enumerate(data)
        ]
    return _parse_tsv(p)


def _parse_tsv(p: Path) -> list[ConsensusFeature]:
    lines = [
        ln for ln in p.read_text(encoding="utf-8").splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    if not lines:
        raise UserFeatureError(f"{p}: empty feature file")
    header = [h.strip().lower() for h in lines[0].split("\t")]
    if "name" not in header or "start" not in header or "end" not in header:
        raise UserFeatureError(
            f"{p}: TSV must have a header line with at least name/start/end "
            f"(got: {header})"
        )
    out = []
    for i, ln in enumerate(lines[1:], start=2):
        vals = ln.split("\t")
        rec = {header[j]: vals[j].strip() for j in range(min(len(header), len(vals)))}
        out.append(feature_from_record(rec, where=f"{p.name}:{i}"))
    return out
