#!/usr/bin/env python3
"""Reproducible demonstration-subset criterion (Application Note).

Criterion (documented, reproducible — NOT a favorable screenshot pick):
the 20 FULL-LENGTH L1HS copies with the highest mean H3K9me3
fold-change over their observed consensus bins, exactly as exported in
``loci.tsv`` (column ``mean_H3K9me3 fold change``). This is an
EXPLORATORY selection chosen to demonstrate the subset-comparison
workflow; it is not presented as a confirmatory biological finding.

Usage:
  python select_demo_loci.py <results_dir> <out_selection.yaml> [N=20]

Reads <results_dir>/loci.tsv + analysis.yaml (provenance) and writes a
tescope-locus-selection/1 YAML importable in the GUI's Loci tab.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

MEAN_COL = "mean_H3K9me3 fold change"


def main() -> None:
    results = Path(sys.argv[1])
    out = Path(sys.argv[2])
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    loci = pd.read_csv(results / "loci.tsv", sep="\t")
    cfg = yaml.safe_load((results / "analysis.yaml").read_text())
    full = loci[loci["completeness"] == "full_length"].copy()
    full = full[full[MEAN_COL].notna()]
    top = full.sort_values(MEAN_COL, ascending=False).head(n)
    doc = {
        "schema": "tescope-locus-selection/1",
        "assembly": cfg["assembly"],
        "repeat": cfg["repeat"]["name"],
        "coordinates": "0-based half-open",
        "annotation_provenance": cfg.get("annotation_provenance", {}),
        "n_eligible_at_export": int(len(loci)),
        "selection_criterion": (
            f"top {n} full_length copies by '{MEAN_COL}' "
            "(mean over observed consensus bins, loci.tsv); exploratory "
            "workflow demonstration, not a confirmatory finding"
        ),
        "locus_ids": top["locus_id"].tolist(),
    }
    out.write_text(yaml.safe_dump(doc, sort_keys=False))
    print(f"wrote {out}: {len(top)} loci "
          f"(mean range {top[MEAN_COL].min():.3g}–{top[MEAN_COL].max():.3g})")


if __name__ == "__main__":
    main()
