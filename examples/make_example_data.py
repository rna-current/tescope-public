#!/usr/bin/env python
"""Generate the TEscope demo dataset (deterministic, seeded).

Creates, in ./example-data/:

- ``demo.fa.out``          RepeatMasker annotation: 60 copies of "L1demo"
                           (6 kb consensus) on two chromosomes, realistic
                           fragment-length structure, both strands
- ``demo.chrom.sizes``
- ``rna_ctrl_rep{1,2}.bedGraph`` / ``rna_treat_rep{1,2}.bedGraph``
                           coverage with a 3'-biased profile; "treatment"
                           adds a bump over a 5' region
- ``h3k9me3.bedGraph``     broad signal enriched on full-length copies
- ``peaks.narrowPeak``     peaks over the 5' region of near-full copies
- ``sites.bed``            point events clustered near consensus 5.5-6 kb
- ``analysis.yaml``        ready-to-run project file

Usage:
    python make_example_data.py [outdir]
    tescope reference add --assembly demo --repeatmasker example-data/demo.fa.out \
        --chrom-sizes example-data/demo.chrom.sizes
    tescope run example-data/analysis.yaml
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

CONS_LEN = 6000
CHROMS = {"chr1": 2_000_000, "chr2": 1_500_000}
SEED = 20260827


def make_copies(rng: np.random.Generator):
    """Return list of (chrom, gs, ge, strand, cons_start, cons_end)."""
    copies = []
    cursor = {c: 50_000 for c in CHROMS}
    for i in range(60):
        chrom = "chr1" if i % 3 else "chr2"
        r = rng.random()
        if r < 0.15:  # full length
            cs, ce = 0, CONS_LEN
        elif r < 0.70:  # 5'-truncated (L1-typical 3' fragments)
            ce = CONS_LEN
            cs = int(rng.integers(1500, 5700))
        elif r < 0.85:  # 3'-truncated
            cs = 0
            ce = int(rng.integers(400, 4500))
        else:  # internal
            cs = int(rng.integers(500, 4000))
            ce = cs + int(rng.integers(300, 1500))
        length = ce - cs  # 1:1 genomic:consensus for simplicity
        gs = cursor[chrom] + int(rng.integers(5_000, 20_000))
        ge = gs + length
        cursor[chrom] = ge
        strand = "+" if rng.random() < 0.5 else "-"
        copies.append((chrom, gs, ge, strand, cs, ce))
    return copies


def write_out(path: Path, copies):
    lines = [
        "   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat",
        "score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID",
        "",
    ]
    for i, (chrom, gs, ge, strand, cs, ce) in enumerate(copies, 1):
        left = CONS_LEN - ce
        qleft = CHROMS[chrom] - ge
        if strand == "+":
            r1, r2, r3, s = str(cs + 1), str(ce), f"({left})", "+"
        else:
            r1, r2, r3, s = f"({left})", str(ce), str(cs + 1), "C"
        lines.append(
            f"{2000 + i:>5}  5.0  1.0  1.0  {chrom:<8} {gs + 1:>8} {ge:>8} "
            f"{'(' + str(qleft) + ')':>12} {s} {'L1demo':<10} {'LINE/L1':<18} "
            f"{r1:>6} {r2:>6} {r3:>7} {i}"
        )
    path.write_text("\n".join(lines) + "\n")


def profile_value(cons_pos: np.ndarray, condition: str, rng) -> np.ndarray:
    """Consensus-position-dependent signal: 3' bias + treatment 5' bump."""
    base = 1.0 + 4.0 * (cons_pos / CONS_LEN) ** 2  # 3'-biased (poly-A proximal)
    if condition == "treat":
        base = base + 6.0 * np.exp(-((cons_pos - 900) ** 2) / (2 * 250**2))
    noise = rng.normal(0, 0.15, size=len(cons_pos))
    return np.maximum(base + noise, 0.0)


def write_bedgraph(path: Path, copies, condition: str, rng, broad_full_only=False):
    with open(path, "w") as fh:
        for chrom, gs, ge, strand, cs, ce in copies:
            L = ge - gs
            i = np.arange(L)
            cons = (cs + i) if strand == "+" else (cs + (L - 1 - i))
            if broad_full_only:
                vals = np.full(L, 3.0 if (ce - cs) > 0.8 * CONS_LEN else 0.5)
                vals = vals + rng.normal(0, 0.1, L)
            else:
                vals = profile_value(cons.astype(float), condition, rng)
            # run-length encode to keep files small
            start = 0
            for k in range(1, L + 1):
                if k == L or abs(vals[k] - vals[start]) > 1e-12:
                    fh.write(f"{chrom}\t{gs + start}\t{gs + k}\t{vals[start]:.3f}\n")
                    start = k


def write_peaks(path: Path, copies, rng):
    with open(path, "w") as fh:
        n = 0
        for chrom, gs, ge, strand, cs, ce in copies:
            if cs > 1200:  # peak sits at consensus ~800-1000; needs 5' end
                continue
            if rng.random() < 0.3:
                continue
            # genomic interval corresponding to consensus [800, 1000)
            if strand == "+":
                a, b = gs + (800 - cs), gs + (1000 - cs)
            else:
                a, b = ge - (1000 - cs), ge - (800 - cs)
            a, b = max(a, gs), min(b, ge)
            if b <= a:
                continue
            n += 1
            fh.write(f"{chrom}\t{a}\t{b}\tpeak{n}\t500\t.\t{8 + rng.random() * 4:.2f}\t5.0\t3.0\t{(b - a) // 2}\n")


def write_sites(path: Path, copies, rng):
    with open(path, "w") as fh:
        n = 0
        for chrom, gs, ge, strand, cs, ce in copies:
            span = ce - cs
            for _ in range(int(rng.integers(0, 4))):
                cons_target = int(rng.normal(5600, 200))
                if not (cs <= cons_target < ce):
                    continue
                if strand == "+":
                    pos = gs + (cons_target - cs)
                else:
                    pos = ge - 1 - (cons_target - cs)
                n += 1
                fh.write(f"{chrom}\t{pos}\t{pos + 1}\tsite{n}\t{rng.random():.3f}\t{strand}\n")
            _ = span
    return


ANALYSIS_YAML = """\
assembly: demo
repeat:
  name: L1demo
filters:
  min_consensus_coverage: 0.1
tracks:
  - file: rna_ctrl_rep1.bedGraph
    name: RNA control r1
    group: control
    replicate: 1
  - file: rna_ctrl_rep2.bedGraph
    name: RNA control r2
    group: control
    replicate: 2
  - file: rna_treat_rep1.bedGraph
    name: RNA treated r1
    group: treated
    replicate: 1
  - file: rna_treat_rep2.bedGraph
    name: RNA treated r2
    group: treated
    replicate: 2
  - file: h3k9me3.bedGraph
    name: H3K9me3
  - file: peaks.narrowPeak
    name: 5' peaks
    kind: interval
  - file: sites.bed
    name: modification sites
    kind: point
aggregation:
  bin_size: 50
  statistic: mean
output:
  dir: demo-results
  figure_formats: [png, svg]
"""


def main() -> None:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "example-data")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    copies = make_copies(rng)
    write_out(out / "demo.fa.out", copies)
    (out / "demo.chrom.sizes").write_text(
        "".join(f"{c}\t{ln}\n" for c, ln in CHROMS.items())
    )
    for cond in ("ctrl", "treat"):
        for rep in (1, 2):
            write_bedgraph(
                out / f"rna_{cond}_rep{rep}.bedGraph", copies, cond,
                np.random.default_rng(SEED + hash((cond, rep)) % 10_000),
            )
    write_bedgraph(out / "h3k9me3.bedGraph", copies, "ctrl",
                   np.random.default_rng(SEED + 77), broad_full_only=True)
    write_peaks(out / "peaks.narrowPeak", copies, np.random.default_rng(SEED + 88))
    write_sites(out / "sites.bed", copies, np.random.default_rng(SEED + 99))
    (out / "analysis.yaml").write_text(ANALYSIS_YAML)
    print(f"demo dataset written to {out}/")
    print("next:")
    print(f"  tescope reference add --assembly demo --repeatmasker {out}/demo.fa.out "
          f"--chrom-sizes {out}/demo.chrom.sizes")
    print(f"  tescope run {out}/analysis.yaml")


if __name__ == "__main__":
    main()
