#!/usr/bin/env python
"""Build the portable TEscope acceptance-test pack (SYNTHETIC data).

Produces ``acceptance-pack/`` — a self-contained, offline, relocatable
directory that exercises the application's core functionality with
deterministic inputs and KNOWN expected outputs:

  annotation (RepeatMasker .out, 60 synthetic "L1demo" copies)
  continuous signal (bigWig + bedGraph), interval peaks (narrowPeak)
  small BAM with controlled NH mapping cases (+ .bai)
  user consensus FASTA + user feature table (sequence/annotation layer)
  three ready-to-run configs (core / mapping-classes / sequence)
  expected/ canonical TSVs + SHA256SUMS  (generated here, from this
      TEscope version)
  verify.py — reruns all three analyses offline and byte-compares
      the canonical tables against expected/
  QUICKSTART.md — CLI + GUI instructions

Everything is synthetic and labeled synthetic; this pack tests
functionality, not biology. Regenerate with:
    .venv/bin/python examples/make_acceptance_pack.py [outdir]
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_example_data as demo  # noqa: E402  (deterministic generators)

SEED = 20261002
CONS_LEN = demo.CONS_LEN


def write_bigwig(path: Path, bedgraph: Path, chrom_sizes: dict[str, int]):
    import pyBigWig

    bw = pyBigWig.open(str(path), "w")
    bw.addHeader(sorted(chrom_sizes.items()))
    entries: dict[str, list] = {}
    for line in bedgraph.read_text().splitlines():
        c, s, e, v = line.split()
        entries.setdefault(c, []).append((int(s), int(e), float(v)))
    for c in sorted(entries):
        rows = sorted(entries[c])
        bw.addEntries([c] * len(rows), [r[0] for r in rows],
                      ends=[r[1] for r in rows], values=[r[2] for r in rows])
    bw.close()


def write_bam(path: Path, copies, chrom_sizes):
    """Controlled NH cases on the first full-length copy:
    40 unique reads (NH=1) across its 5' kb, 30 multi reads (NH=2) on
    its body, 10 unique reads on an intergenic position (outside any
    copy)."""
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": ln} for c, ln in
                     sorted(chrom_sizes.items())]}
    ref_ids = {c: i for i, (c, _ln) in enumerate(sorted(chrom_sizes.items()))}
    full = next(c for c in copies if c[4] <= 5 and CONS_LEN - c[5] <= 5)
    chrom, gs, _ge, _strand, _cs, _ce = full
    reads = []
    for i in range(40):  # unique over the copy 5' region
        reads.append((f"u{i}", gs + 20 + i * 20, 1))
    for i in range(30):  # multi on the body
        reads.append((f"m{i}", gs + 2000 + i * 30, 2))
    for i in range(10):  # unique, intergenic
        reads.append((f"bg{i}", 10_000 + i * 50, 1))
    tmp = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(tmp, "wb", header=header) as bf:
        for name, pos, nh in sorted(reads, key=lambda r: r[1]):
            a = pysam.AlignedSegment(bf.header)
            a.query_name = name
            a.query_sequence = "A" * 50
            a.flag = 0
            a.reference_id = ref_ids[chrom]
            a.reference_start = pos
            a.mapping_quality = 255 if nh == 1 else 1
            a.cigarstring = "50M"
            a.set_tag("NH", nh)
            bf.write(a)
    pysam.sort("-o", str(path), str(tmp))
    tmp.unlink()
    pysam.index(str(path))


def write_fasta_and_features(fa: Path, feats: Path):
    rng = np.random.default_rng(SEED)
    seq = "".join(rng.choice(list("ACGT"), size=CONS_LEN))
    fa.write_text(">L1demo#LINE/L1 synthetic demo consensus\n" + "\n".join(
        seq[i:i + 80] for i in range(0, CONS_LEN, 80)) + "\n")
    feats.write_text(
        "name\tstart\tend\tlabel\n"
        "FiveUTR\t1\t900\t5'UTR (synthetic)\n"
        "ORFdemo\t901\t5400\tORF (synthetic)\n"
        "ThreeUTR\t5401\t6000\t3'UTR (synthetic)\n")


CONFIG_CORE = """# acceptance: core profiles (continuous + interval tracks)
tescope_version: {ver}
assembly: demo
repeat:
  name: L1demo
tracks:
- file: inputs/signal.bw
  name: signal bigWig
- file: inputs/h3k9me3.bedGraph
  name: broad bedGraph
- file: inputs/peaks.narrowPeak
  name: demo peaks
aggregation:
  bin_size: 100
output:
  dir: out/core
  table_formats: [tsv]
  figure_formats: [png]
"""

CONFIG_MAPPING = """# acceptance: BAM NH mapping classes (unique/multi gating)
tescope_version: {ver}
assembly: demo
repeat:
  name: L1demo
tracks:
- file: inputs/reads.bam
  name: reads all
- file: inputs/reads.bam
  name: reads unique
  mapping_class: unique
  mapping_method: nh
- file: inputs/reads.bam
  name: reads multi
  mapping_class: multi
  mapping_method: nh
aggregation:
  bin_size: 200
output:
  dir: out/mapping
  table_formats: [tsv]
  figure_formats: [png]
"""

CONFIG_SEQ = """# acceptance: user consensus sequence + user features
tescope_version: {ver}
assembly: demo
repeat:
  name: L1demo
tracks:
- file: inputs/signal.bw
  name: signal bigWig
aggregation:
  bin_size: 100
annotations:
  enabled: true
  user_features_file: inputs/user_features.tsv
sequence:
  enabled: true
  source: fasta
  path: inputs/L1demo_consensus.fa
output:
  dir: out/sequence
  table_formats: [tsv]
  figure_formats: [png]
"""

VERIFY = '''#!/usr/bin/env python
"""Offline acceptance verification: reruns the three analyses and
byte-compares every canonical TSV against expected/. Requires only an
installed `tescope` CLI; run from the pack directory:
    python verify.py
"""
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

PACK = Path(__file__).resolve().parent
CACHE = PACK / ".tescope-cache"


def sh(*cmd):
    print("+", " ".join(map(str, cmd)))
    subprocess.run(list(map(str, cmd)), check=True, cwd=PACK)


def main():
    tescope = shutil.which("tescope") or sys.exit("tescope CLI not on PATH")
    shutil.rmtree(PACK / "out", ignore_errors=True)
    shutil.rmtree(CACHE, ignore_errors=True)
    env_cache = ["--cache-dir", str(CACHE)]
    sh(tescope, "reference", "add", "--assembly", "demo",
       "--repeatmasker", "inputs/demo.fa.out",
       "--chrom-sizes", "inputs/demo.chrom.sizes", *env_cache)
    for cfg in ("config_core.yaml", "config_mapping.yaml",
                "config_sequence.yaml"):
        sh(tescope, "run", cfg, *env_cache)
    bad = 0
    for sub in ("core", "mapping", "sequence"):
        if not (PACK / "out" / sub / "qc.tsv").exists():
            print("MISSING:", sub + "/qc.tsv"); bad += 1
    sums = (PACK / "expected/SHA256SUMS").read_text().splitlines()
    for line in sums:
        digest, rel = line.split(maxsplit=1)
        got = PACK / "out" / rel
        if not got.exists():
            print("MISSING:", rel); bad += 1; continue
        d = hashlib.sha256(got.read_bytes()).hexdigest()
        if d != digest:
            print("MISMATCH:", rel); bad += 1
    n = len(sums)
    print(f"{n - bad}/{n} canonical tables match expected checksums")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
'''

QUICKSTART = """# TEscope acceptance pack — quick start (SYNTHETIC data)

Everything here is deterministic synthetic test data ("L1demo"); it
verifies that TEscope works, not any biological claim. The pack is
relocatable and runs fully offline.

## CLI (installed wheel: `pip install tescope` or `uv tool install …`)

    cd <this directory>
    python verify.py

This registers the bundled synthetic reference, runs the three
configurations (core profiles, BAM NH mapping classes, user
sequence/features), and byte-compares every canonical TSV against
`expected/` (exit 0 = pass).

## GUI (installed desktop app)

1. Reference genome ▸ Import custom… → `inputs/demo.fa.out` +
   `inputs/demo.chrom.sizes` (assembly name: demo).
2. Select repeat `L1demo`; add `inputs/signal.bw`,
   `inputs/h3k9me3.bedGraph`, `inputs/peaks.narrowPeak`; Run analysis.
3. Expect: 60 copies; a 3'-biased continuous profile; peak occurrence
   concentrated in the 5' region of near-full-length copies.
4. Mapping classes: load `config_mapping.yaml` (File ▸ Open config) —
   unique (NH=1) coverage sits on the 5' kb of one full-length copy,
   multi (NH>1) on its body.
5. Sequence/features: open `config_sequence.yaml` — the user consensus
   (6,000 bp) displays with three synthetic feature intervals.
6. Export figure/data; File ▸ Save config; About ▸ Help for offline help.

Expected numerical outputs (canonical TSVs) and their SHA-256 sums are
in `expected/`; figures are style-dependent and not checksummed.
"""


def main() -> None:
    # resolve(): expected-output subprocesses run with cwd=out, so a
    # relative out dir would mis-place --cache-dir inside the pack
    out = (Path(sys.argv[1]) if len(sys.argv) > 1
           else HERE / "acceptance-pack").resolve()
    if out.exists():
        shutil.rmtree(out)
    inputs = out / "inputs"
    inputs.mkdir(parents=True)

    rng = np.random.default_rng(demo.SEED)
    copies = demo.make_copies(rng)
    demo.write_out(inputs / "demo.fa.out", copies)
    (inputs / "demo.chrom.sizes").write_text(
        "".join(f"{c}\t{ln}\n" for c, ln in demo.CHROMS.items()))
    demo.write_bedgraph(inputs / "h3k9me3.bedGraph", copies, "h3k9me3",
                        rng, broad_full_only=True)
    demo.write_bedgraph(inputs / "signal.bedGraph", copies, "ctrl", rng)
    demo.write_peaks(inputs / "peaks.narrowPeak", copies, rng)
    write_bigwig(inputs / "signal.bw", inputs / "signal.bedGraph",
                 demo.CHROMS)
    (inputs / "signal.bedGraph").unlink()  # bigWig is the shipped form
    write_bam(inputs / "reads.bam", copies, demo.CHROMS)
    write_fasta_and_features(inputs / "L1demo_consensus.fa",
                             inputs / "user_features.tsv")

    from tescope import __version__
    for name, tmpl in (("config_core.yaml", CONFIG_CORE),
                       ("config_mapping.yaml", CONFIG_MAPPING),
                       ("config_sequence.yaml", CONFIG_SEQ)):
        (out / name).write_text(tmpl.format(ver=__version__))
    (out / "verify.py").write_text(VERIFY)
    (out / "QUICKSTART.md").write_text(QUICKSTART)

    # ---- generate expected outputs with THIS TEscope version ---------
    tescope = HERE.parent / ".venv/bin/tescope"
    cache = out / ".build-cache"
    subprocess.run([str(tescope), "reference", "add", "--assembly", "demo",
                    "--repeatmasker", "inputs/demo.fa.out",
                    "--chrom-sizes", "inputs/demo.chrom.sizes",
                    "--cache-dir", str(cache)], check=True, cwd=out)
    for cfg in ("config_core.yaml", "config_mapping.yaml",
                "config_sequence.yaml"):
        subprocess.run([str(tescope), "run", cfg, "--cache-dir", str(cache)],
                       check=True, cwd=out)
    exp = out / "expected"
    exp.mkdir()
    lines = []
    for tsv in sorted((out / "out").rglob("*.tsv")):
        rel = tsv.relative_to(out / "out")
        if rel.name == "qc.tsv":
            # qc embeds run timestamps + absolute input paths
            # (provenance, environment-specific) — verified for
            # existence, not byte identity
            continue
        dest = exp / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tsv, dest)
        digest = hashlib.sha256(tsv.read_bytes()).hexdigest()
        lines.append(f"{digest}  {rel}")
    (exp / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    shutil.rmtree(out / "out")
    shutil.rmtree(cache)
    print(f"acceptance pack built: {out} "
          f"({len(lines)} expected tables; "
          f"{sum(f.stat().st_size for f in out.rglob('*') if f.is_file())/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
