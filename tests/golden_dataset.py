"""Deterministic synthetic golden dataset.

One artificial chromosome ``chrT`` (100 kb) carrying copies of an invented
repeat ``TST1`` (consensus length 1000) plus one ``TST2`` copy:

====== ================== ======== ============== =========================
copy   genomic (0-based)  strand   consensus      completeness (tol 50 bp)
====== ================== ======== ============== =========================
F      10000-11000        +        [0, 1000)      full_length
R      20000-21000        -        [0, 1000)      full_length
P5     30000-30400        +        [0, 400)       three_prime_truncated
P3R    40000-40300        -        [700, 1000)    five_prime_truncated
MID    50000-50200        +        [400, 600)     internal_fragment
SHORT  60000-60050        +        [100, 150)     internal_fragment
STAR   90000-90100        +        [0, 100)       secondary hit ('*')
====== ================== ======== ============== =========================

The continuous signal is constructed so that, under **correct** strand-aware
consensus mapping, every genomic base carries a value equal to its consensus
position. Aggregated with bin size B, every covered bin must therefore have
mean == median == bin center (b*B + (B-1)/2), for every copy composition —
and any orientation or offset bug shifts or mirrors the ramp.

Peaks:
- ``peak_F``  chrT:10250-10300  -> consensus [250, 300)   (forward copy)
- ``peak_R``  chrT:20100-20200  -> consensus [800, 900)   (reverse copy!)
- ``peak_out`` chrT:70000-70100 -> projects nowhere

Point sites:
- chrT:10050 -> consensus 50.5 (copy F), chrT:20950 -> consensus 49.5 (copy R)
"""

from __future__ import annotations

import os
from pathlib import Path

CHROM = "chrT"
CHROM_LEN = 100_000
CONS_LEN = 1000

# name -> (genome_start, genome_end, strand, cons_start, cons_end, secondary)
COPIES = {
    "F": (10_000, 11_000, "+", 0, 1000, False),
    "R": (20_000, 21_000, "-", 0, 1000, False),
    "P5": (30_000, 30_400, "+", 0, 400, False),
    "P3R": (40_000, 40_300, "-", 700, 1000, False),
    "MID": (50_000, 50_200, "+", 400, 600, False),
    "SHORT": (60_000, 60_050, "+", 100, 150, False),
    "STAR": (90_000, 90_100, "+", 0, 100, True),
}

EXPECTED_COMPLETENESS = {
    "F": "full_length",
    "R": "full_length",
    "P5": "three_prime_truncated",
    "P3R": "five_prime_truncated",
    "MID": "internal_fragment",
    "SHORT": "internal_fragment",
    "STAR": "three_prime_truncated",
}


def per_base_values() -> list[tuple[int, float]]:
    """(genomic position, value) for every base inside a TST1 copy.

    value == the consensus position of the base under correct mapping.
    """
    out = []
    for gs, ge, strand, cs, ce, _sec in COPIES.values():
        L = ge - gs
        for i in range(L):
            cons = (cs + i) if strand == "+" else (cs + (L - 1 - i))
            out.append((gs + i, float(cons)))
    out.sort()
    return out


def _out_line(sw, div, chrom, g1, g2, strand, name, cf, r1, r2, r3, rm_id, star=False):
    tail = f" {rm_id}" + (" *" if star else "")
    return (
        f"{sw:>5} {div:4.1f}  0.0  0.0  {chrom:<8} {g1:>8} {g2:>8} "
        f"{'(' + str(CHROM_LEN - g2) + ')':>10} {strand} {name:<10} {cf:<18} {r1:>6} {r2:>6} {r3:>7}{tail}"
    )


def write_golden_out(path: str | os.PathLike) -> None:
    lines = [
        "   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat",
        "score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID",
        "",
    ]
    rm_id = 0
    for cname, (gs, ge, strand, cs, ce, sec) in COPIES.items():
        rm_id += 1
        g1, g2 = gs + 1, ge  # 1-based inclusive
        left = CONS_LEN - ce
        if strand == "+":
            r1, r2, r3 = str(cs + 1), str(ce), f"({left})"
            s = "+"
        else:
            r1, r2, r3 = f"({left})", str(ce), str(cs + 1)
            s = "C"
        lines.append(
            _out_line(1000 + rm_id, 1.0, CHROM, g1, g2, s, "TST1", "LINE/TST", r1, r2, r3, rm_id, star=sec)
        )
    # one TST2 copy, different class/family
    lines.append(
        _out_line(500, 5.0, CHROM, 80_001, 80_500, "+", "TST2", "SINE/TST2", "1", "500", "(0)", 99)
    )
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_chrom_sizes(path: str | os.PathLike) -> None:
    Path(path).write_text(f"{CHROM}\t{CHROM_LEN}\n", encoding="utf-8")


def write_bedgraph(path: str | os.PathLike, chrom: str = CHROM) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for pos, val in per_base_values():
            fh.write(f"{chrom}\t{pos}\t{pos + 1}\t{val}\n")


def write_bigwig(path: str | os.PathLike) -> str:
    """Write the signal as bigWig with whichever backend is available.

    Returns the backend used ("pybigwig" or "pybigtools").
    """
    entries = per_base_values()
    starts = [p for p, _ in entries]
    values = [v for _, v in entries]
    try:
        import pyBigWig

        bw = pyBigWig.open(os.fspath(path), "w")
        bw.addHeader([(CHROM, CHROM_LEN)])
        bw.addEntries(
            [CHROM] * len(starts),
            starts,
            ends=[s + 1 for s in starts],
            values=values,
        )
        bw.close()
        return "pybigwig"
    except ImportError:
        pass
    import pybigtools

    # pybigtools' Rust URL detection misreads Windows drive-letter paths
    # ("C:\...") as URL schemes and refuses to write; open by bare filename
    # from the target directory instead (no colon -> no scheme).
    p = Path(path)
    cwd = os.getcwd()
    try:
        os.chdir(p.parent)
        bw = pybigtools.open(p.name, "w")
        bw.write(
            {CHROM: CHROM_LEN},
            ((CHROM, s, s + 1, float(v)) for s, v in zip(starts, values)),
        )
    finally:
        os.chdir(cwd)
    return "pybigtools"


PEAKS = {
    "peak_F": (10_250, 10_300),
    "peak_R": (20_100, 20_200),
    "peak_out": (70_000, 70_100),
}
# expected consensus projections: {peak_name: (cons_start, cons_end)}
EXPECTED_PEAK_PROJECTIONS = {"peak_F": (250, 300), "peak_R": (800, 900)}

SITES = [10_050, 20_950]
EXPECTED_SITE_CONS = [50.5, 49.5]


def write_peaks_bed(path: str | os.PathLike, chrom: str = CHROM) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, (s, e) in PEAKS.items():
            fh.write(f"{chrom}\t{s}\t{e}\t{name}\t100\t+\n")


def write_narrowpeak(path: str | os.PathLike, chrom: str = CHROM) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, (s, e) in PEAKS.items():
            fh.write(f"{chrom}\t{s}\t{e}\t{name}\t100\t.\t7.5\t3.0\t2.0\t{(e - s) // 2}\n")


def write_sites_bed(path: str | os.PathLike, chrom: str = CHROM) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for i, pos in enumerate(SITES):
            fh.write(f"{chrom}\t{pos}\t{pos + 1}\tsite{i}\t{float(i + 1)}\t+\n")


def build_all(directory: str | os.PathLike) -> dict[str, str]:
    """Materialize the full golden dataset into a directory."""
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    paths = {
        "out": d / "golden.out",
        "chrom_sizes": d / "chrom.sizes",
        "bedgraph": d / "signal.bedGraph",
        "bedgraph_alias": d / "signal_alias.bedGraph",
        "bigwig": d / "signal.bw",
        "peaks_bed": d / "peaks.bed",
        "narrowpeak": d / "peaks.narrowPeak",
        "sites": d / "sites.bed",
    }
    write_golden_out(paths["out"])
    write_chrom_sizes(paths["chrom_sizes"])
    write_bedgraph(paths["bedgraph"])
    write_bedgraph(paths["bedgraph_alias"], chrom="T")  # tests chr-name aliasing
    write_bigwig(paths["bigwig"])
    write_peaks_bed(paths["peaks_bed"])
    write_narrowpeak(paths["narrowpeak"])
    write_sites_bed(paths["sites"])
    return {k: str(v) for k, v in paths.items()}
