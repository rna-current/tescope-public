#!/usr/bin/env python
"""SYNTHETIC NANOPORE DEMONSTRATION generator — NOT BIOLOGICAL EVIDENCE.

Every file and figure produced from this dataset must carry the label
"SYNTHETIC NANOPORE DEMONSTRATION". Nothing here is an experimental
finding: the data are generated with a fixed seed and a KNOWN injected
ground truth, so the demo shows that TEscope's consensus projection and
mapping-classification machinery recover exactly what was put in — a
software validation example.

Design (chromosome chrN, 400 kb; repeat "L1nano", consensus 6,000 bp):
- 20 full-length copies (consensus 0-6000) and 20 5'-truncated copies
  (consensus 2500-6000), both strands.
- Long-read alignment BAM: reads span large parts of each copy;
  NH tags written explicitly — reads on full-length copies NH=1
  ("unique-like"), reads on truncated copies NH=4 with MAPQ 0
  ("multi-like"). This is a SIMULATED mapping structure.
- Modification-probability bedGraph: baseline probability 0.05 with a
  KNOWN injected high-probability region at consensus 4000-4500
  (p ~ 0.85) present on every copy that covers it, plus seeded noise.
- GROUND TRUTH: the recovered consensus profile must show the elevated
  region at 4000-4500 and coverage restricted to each copy's true
  consensus footprint (nothing 5' of 2500 for truncated copies).

Usage:
    python examples/make_nanopore_demo.py [outdir]     (default: nanopore-demo)
    tescope reference add --assembly nanodemo \
        --repeatmasker nanopore-demo/nanodemo.fa.out \
        --chrom-sizes nanopore-demo/nanodemo.chrom.sizes
    tescope run nanopore-demo/analysis_modprob.yaml
    tescope run nanopore-demo/analysis_mapping.yaml   (needs tescope[bam])
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SEED = 20260904
CHROM, CHROM_LEN = "chrN", 400_000
CONS_LEN = 6000
MOD_REGION = (4000, 4500)  # injected ground truth (consensus coordinates)
LABEL = "SYNTHETIC NANOPORE DEMONSTRATION - NOT BIOLOGICAL EVIDENCE"


def make_copies(rng):
    copies, cursor = [], 10_000
    for i in range(40):
        full = i < 20
        cs, ce = (0, CONS_LEN) if full else (2500, CONS_LEN)
        length = ce - cs
        gs = cursor + int(rng.integers(500, 2_000))
        cursor = gs + length
        strand = "+" if i % 2 == 0 else "-"
        copies.append((gs, gs + length, strand, cs, ce, full))
    return copies


def write_out(path, copies):
    lines = [
        "   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat",
        "score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID",
        "",
    ]
    for i, (gs, ge, strand, cs, ce, _f) in enumerate(copies, 1):
        left = CONS_LEN - ce
        qleft = CHROM_LEN - ge
        if strand == "+":
            r1, r2, r3, sflag = str(cs + 1), str(ce), f"({left})", "+"
        else:
            r1, r2, r3, sflag = f"({left})", str(ce), str(cs + 1), "C"
        lines.append(
            f"{3000 + i:>5}  2.0  0.5  0.5  {CHROM:<8} {gs + 1:>8} {ge:>8} "
            f"{'(' + str(qleft) + ')':>12} {sflag} {'L1nano':<10} {'LINE/L1':<18} "
            f"{r1:>6} {r2:>6} {r3:>7} {i}"
        )
    Path(path).write_text("\n".join(lines) + "\n")


def write_modprob_bedgraph(path, copies, rng):
    """Per-base modification probability with the injected region."""
    with open(path, "w") as fh:
        for gs, ge, strand, cs, ce, _f in copies:
            L = ge - gs
            i = np.arange(L)
            cons = (cs + i) if strand == "+" else (cs + (L - 1 - i))
            p = np.full(L, 0.05) + rng.normal(0, 0.01, L)
            inj = (cons >= MOD_REGION[0]) & (cons < MOD_REGION[1])
            p[inj] = 0.85 + rng.normal(0, 0.03, int(inj.sum()))
            p = np.clip(p, 0.0, 1.0)
            start = 0
            for k in range(1, L + 1):
                if k == L or abs(p[k] - p[start]) > 0.02:
                    fh.write(f"{CHROM}\t{gs + start}\t{gs + k}\t{p[start]:.3f}\n")
                    start = k


def write_longread_bam(path, copies, rng):
    """Long reads with explicit SIMULATED NH tags (unique vs multi-like)."""
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    reads = []
    for ci, (gs, ge, strand, _cs, _ce, full) in enumerate(copies):
        n_reads = int(rng.integers(6, 12))
        for ri in range(n_reads):
            span = int(rng.integers(1_500, min(4_000, ge - gs)))
            start = int(rng.integers(gs, max(ge - span, gs + 1)))
            nh = 1 if full else 4
            mapq = 60 if full else 0
            reads.append((start, span, mapq, nh, strand == "-", f"r{ci}_{ri}"))
    reads.sort()
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for start, span, mapq, nh, rev, name in reads:
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = "A" * span
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = mapq
            a.cigarstring = f"{span}M"
            a.flag = 0x10 if rev else 0
            a.set_tag("NH", nh)
            bam.write(a)
    pysam.index(str(path))


MODPROB_YAML = """\
# {label}
assembly: nanodemo
repeat: {{name: L1nano}}
tracks:
  - {{file: modprob.bedGraph, name: mod probability (synthetic)}}
aggregation: {{bin_size: 50, statistic: mean}}
style:
  preset: publication
  representation: filled
  band: iqr
  title: "{label}: modification probability on L1nano consensus"
output: {{dir: nanodemo-modprob, figure_formats: [png, svg]}}
"""

MAPPING_YAML = """\
# {label}
assembly: nanodemo
repeat: {{name: L1nano}}
tracks:
  - {{file: longreads.bam, name: unique reads (NH=1), mapping_class: unique, mapping_method: nh}}
  - {{file: longreads.bam, name: multi-like reads (NH>1), mapping_class: multi, mapping_method: nh}}
aggregation: {{bin_size: 50, statistic: mean}}
style:
  preset: publication
  series:
    unique reads (NH=1): {{color: "#0072B2"}}
    multi-like reads (NH>1): {{color: "#D55E00", line_style: dashed}}
  title: "{label}: unique vs multi-like long-read coverage"
output: {{dir: nanodemo-mapping, figure_formats: [png, svg]}}
"""


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "nanopore-demo")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    copies = make_copies(rng)
    write_out(out / "nanodemo.fa.out", copies)
    (out / "nanodemo.chrom.sizes").write_text(f"{CHROM}\t{CHROM_LEN}\n")
    write_modprob_bedgraph(out / "modprob.bedGraph", copies, rng)
    try:
        write_longread_bam(out / "longreads.bam", copies, rng)
        bam_ok = True
    except ImportError:
        bam_ok = False
    (out / "analysis_modprob.yaml").write_text(MODPROB_YAML.format(label=LABEL))
    (out / "analysis_mapping.yaml").write_text(MAPPING_YAML.format(label=LABEL))
    (out / "GROUND_TRUTH.md").write_text(
        f"# {LABEL}\n\n"
        f"- seed: {SEED}\n"
        f"- injected modification region: consensus {MOD_REGION[0]}-{MOD_REGION[1]} (p~0.85 vs 0.05 baseline)\n"
        f"- copies: 20 full-length (0-{CONS_LEN}) + 20 5'-truncated (2500-{CONS_LEN})\n"
        f"- SIMULATED mapping: full-length reads NH=1/MAPQ60; truncated-copy reads NH=4/MAPQ0\n"
        f"- expected recovery: elevated profile only in {MOD_REGION[0]}-{MOD_REGION[1]}; no coverage 5' of 2500 for truncated copies;\n"
        f"  unique-vs-multi overlay separates the two simulated read classes exactly\n")
    print(f"{LABEL}\nwritten to {out}/ (BAM: {'yes' if bam_ok else 'skipped - pysam unavailable'})")


if __name__ == "__main__":
    main()
