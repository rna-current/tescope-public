#!/usr/bin/env python
"""Regenerate tests/golden/expected_profile.tsv (the regression snapshot).

Run this ONLY when a calculation change is intentional, and record why in
the commit message. CI compares analysis output against this file on all
platforms (Test I: cross-OS numerical reproducibility).
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tests"))

import golden_dataset  # noqa: E402

from tescope.config.model import (  # noqa: E402
    AggregationConfig,
    AnalysisConfig,
    RepeatSelection,
    TrackConfig,
)
from tescope.core.analysis import run_analysis  # noqa: E402
from tescope.core.reference.store import ReferenceStore  # noqa: E402


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        paths = golden_dataset.build_all(Path(tmp) / "golden")
        store = ReferenceStore(cache_dir=Path(tmp) / "cache")
        store.add_custom("testasm", paths["out"], chrom_sizes_path=paths["chrom_sizes"])
        cfg = AnalysisConfig(
            assembly="testasm",
            repeat=RepeatSelection(name="TST1"),
            tracks=[
                TrackConfig(file=paths["bigwig"], name="signal_bw", group="ctrl", replicate=1),
            ],
            aggregation=AggregationConfig(bin_size=100),
        )
        result = run_analysis(cfg, store=store)
        prof = result.continuous_tracks()[0].profile
        out = REPO / "tests" / "golden" / "expected_profile.tsv"
        prof.to_csv(out, sep="\t", index=False)
        print(f"wrote {out} ({len(prof)} bins)")


if __name__ == "__main__":
    main()
