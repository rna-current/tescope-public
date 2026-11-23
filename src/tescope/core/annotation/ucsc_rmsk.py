"""Parser for UCSC ``rmsk`` table dumps (``rmsk.txt.gz``).

17 tab-separated columns:

```
bin swScore milliDiv milliDel milliIns genoName genoStart genoEnd genoLeft
strand repName repClass repFamily repStart repEnd repLeft id
```

Coordinate semantics (UCSC schema):

- ``genoStart``/``genoEnd`` are already 0-based half-open.
- ``repEnd`` is the hit's end on the consensus (1-based inclusive).
- For ``+`` strand hits: ``repStart`` is the consensus begin (1-based) and
  ``repLeft`` is minus the number of consensus bases remaining.
- For ``-`` strand hits: ``repLeft`` is the consensus begin (1-based) and
  ``repStart`` is minus the number of consensus bases remaining.

Percent divergence/deletion/insertion are stored per-mille (``milliDiv``).
"""

from __future__ import annotations

import os

import pandas as pd

from .schema import finalize_frame

_COLUMNS = [
    "bin",
    "swScore",
    "milliDiv",
    "milliDel",
    "milliIns",
    "genoName",
    "genoStart",
    "genoEnd",
    "genoLeft",
    "strand",
    "repName",
    "repClass",
    "repFamily",
    "repStart",
    "repEnd",
    "repLeft",
    "id",
]


def parse_ucsc_rmsk(path: str | os.PathLike) -> pd.DataFrame:
    """Parse a UCSC rmsk table dump (optionally gzipped) to the canonical frame."""
    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=_COLUMNS,
        dtype={
            "swScore": "int64",
            "milliDiv": "int64",
            "milliDel": "int64",
            "milliIns": "int64",
            "genoName": "str",
            "genoStart": "int64",
            "genoEnd": "int64",
            "strand": "str",
            "repName": "str",
            "repClass": "str",
            "repFamily": "str",
            "repStart": "int64",
            "repEnd": "int64",
            "repLeft": "int64",
            "id": "int64",
        },
        usecols=[c for c in _COLUMNS if c not in ("bin", "genoLeft")],
        compression="infer",
    )
    plus = raw["strand"] == "+"
    cons_begin_1 = raw["repStart"].where(plus, raw["repLeft"])
    cons_left = (-raw["repLeft"]).where(plus, -raw["repStart"])

    df = pd.DataFrame(
        {
            "chrom": raw["genoName"],
            "genome_start": raw["genoStart"],
            "genome_end": raw["genoEnd"],
            "strand": raw["strand"],
            "repeat_name": raw["repName"],
            "repeat_class": raw["repClass"],
            "repeat_family": raw["repFamily"],
            "sw_score": raw["swScore"].astype("float64"),
            "pct_div": raw["milliDiv"] / 10.0,
            "pct_del": raw["milliDel"] / 10.0,
            "pct_ins": raw["milliIns"] / 10.0,
            "cons_start": cons_begin_1 - 1,
            "cons_end": raw["repEnd"],
            "cons_left": cons_left,
            "rm_id": raw["id"],
            "is_secondary": False,
        }
    )
    # Defensive: a handful of rmsk rows historically carry inconsistent
    # repeat coordinates (repEnd < begin). Drop them rather than fail.
    ok = (df["cons_end"] >= df["cons_start"] + 0) & (df["cons_start"] >= 0)
    dropped = int((~ok).sum())
    df = df[ok]
    out = finalize_frame(df)
    out.attrs["n_malformed_lines"] = dropped
    return out
