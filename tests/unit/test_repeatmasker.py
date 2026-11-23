"""RepeatMasker parsing: native .out and UCSC rmsk, coordinate semantics."""

from __future__ import annotations

import gzip

import golden_dataset as gd
import pytest

from tescope.core.annotation.repeatmasker import (
    RepeatMaskerParseError,
    parse_any,
    parse_repeatmasker_out,
    sniff_format,
)
from tescope.core.annotation.schema import consensus_lengths, repeat_hierarchy
from tescope.core.annotation.ucsc_rmsk import parse_ucsc_rmsk


def _row(df, gs):
    sub = df[df["genome_start"] == gs]
    assert len(sub) == 1
    return sub.iloc[0]


class TestOutParser:
    def test_parses_all_hits(self, golden_dir):
        df = parse_repeatmasker_out(golden_dir["out"])
        assert len(df) == len(gd.COPIES) + 1  # + one TST2 hit

    def test_coordinates_are_zero_based_half_open(self, golden_dir):
        df = parse_repeatmasker_out(golden_dir["out"])
        for name, (gs, ge, strand, cs, ce, sec) in gd.COPIES.items():
            r = _row(df, gs)
            assert (r["genome_start"], r["genome_end"]) == (gs, ge), name
            assert r["strand"] == strand, name
            assert (r["cons_start"], r["cons_end"]) == (cs, ce), name
            assert r["cons_left"] == gd.CONS_LEN - ce, name
            assert bool(r["is_secondary"]) == sec, name

    def test_reverse_strand_consensus_fields(self, golden_dir):
        # C-strand hits print repeat coords as "(left) end begin"
        df = parse_repeatmasker_out(golden_dir["out"])
        r = _row(df, 40_000)  # P3R
        assert r["strand"] == "-"
        assert r["cons_start"] == 700
        assert r["cons_end"] == 1000
        assert r["cons_left"] == 0

    def test_class_family_split(self, golden_dir):
        df = parse_repeatmasker_out(golden_dir["out"])
        r = _row(df, 10_000)
        assert r["repeat_class"] == "LINE"
        assert r["repeat_family"] == "TST"

    def test_scores_and_divergence_kept(self, golden_dir):
        df = parse_repeatmasker_out(golden_dir["out"])
        r = _row(df, 10_000)
        assert r["sw_score"] > 0
        assert r["pct_div"] == pytest.approx(1.0)

    def test_gz_input(self, golden_dir, tmp_path):
        gz = tmp_path / "golden.out.gz"
        with open(golden_dir["out"], "rb") as fin, gzip.open(gz, "wb") as fout:
            fout.write(fin.read())
        df = parse_repeatmasker_out(gz)
        assert len(df) == len(gd.COPIES) + 1

    def test_rejects_garbage(self, tmp_path):
        p = tmp_path / "bad.out"
        p.write_text("this is not repeatmasker output\n")
        with pytest.raises(RepeatMaskerParseError):
            parse_repeatmasker_out(p)

    def test_consensus_lengths(self, golden_dir):
        df = parse_repeatmasker_out(golden_dir["out"])
        lens = consensus_lengths(df)
        assert lens["TST1"] == gd.CONS_LEN
        assert lens["TST2"] == 500

    def test_hierarchy(self, golden_dir):
        df = parse_repeatmasker_out(golden_dir["out"])
        h = repeat_hierarchy(df)
        tst1 = h[h["repeat_name"] == "TST1"].iloc[0]
        assert tst1["repeat_class"] == "LINE"
        assert tst1["n_copies"] == len(gd.COPIES)


class TestUcscRmsk:
    def _make_rmsk(self, tmp_path):
        """Same copies encoded in UCSC rmsk semantics."""
        lines = []
        rm_id = 0
        for gs, ge, strand, cs, ce, _sec in gd.COPIES.values():
            rm_id += 1
            left = gd.CONS_LEN - ce
            if strand == "+":
                rep_start, rep_end, rep_left = cs + 1, ce, -left
            else:
                rep_start, rep_end, rep_left = -left, ce, cs + 1
            lines.append(
                "\t".join(
                    str(x)
                    for x in [
                        585, 1000, 10, 0, 0, gd.CHROM, gs, ge, -(gd.CHROM_LEN - ge),
                        strand, "TST1", "LINE", "TST", rep_start, rep_end, rep_left, rm_id,
                    ]
                )
            )
        p = tmp_path / "rmsk.txt"
        p.write_text("\n".join(lines) + "\n")
        return p

    def test_matches_out_semantics(self, tmp_path, golden_dir):
        rmsk = parse_ucsc_rmsk(self._make_rmsk(tmp_path))
        out = parse_repeatmasker_out(golden_dir["out"])
        out = out[out["repeat_name"] == "TST1"]
        for gs, _ge, _strand, cs, ce, _sec in gd.COPIES.values():
            a = _row(rmsk, gs)
            b = _row(out, gs)
            for col in ("genome_end", "strand", "cons_start", "cons_end", "cons_left"):
                assert a[col] == b[col], (gs, col)
        assert rmsk["pct_div"].iloc[0] == pytest.approx(1.0)  # milliDiv 10 -> 1.0%

    def test_sniff(self, tmp_path, golden_dir):
        assert sniff_format(golden_dir["out"]) == "repeatmasker-out"
        assert sniff_format(self._make_rmsk(tmp_path)) == "ucsc-rmsk"

    def test_parse_any(self, tmp_path):
        df, fmt = parse_any(self._make_rmsk(tmp_path))
        assert fmt == "ucsc-rmsk"
        assert len(df) == len(gd.COPIES)


class TestRealWorldArtifacts:
    def test_artifact_lines_skipped_and_counted(self, golden_dir, tmp_path):
        # real hg38 .out contains rare hits with consensus begin == end + 1;
        # they must be skipped (counted), not abort the whole parse
        src = open(golden_dir["out"]).read()
        bad = (
            " 1891  11.8  2.8  4.5  chr1      4639264 4639274 (244317148) C "
            " AluSz6         SINE/Alu             (176)  136    137   5612\n"
        )
        p = tmp_path / "with_artifact.out"
        p.write_text(src + bad)
        df = parse_repeatmasker_out(p)
        assert df.attrs["n_malformed_lines"] == 1
        assert len(df) == len(gd.COPIES) + 1

    def test_mostly_garbage_still_fails(self, tmp_path):
        lines = "\n".join(f"garbage line {i} x y z" for i in range(50))
        p = tmp_path / "garbage.out"
        p.write_text(lines + "\n")
        with pytest.raises(RepeatMaskerParseError):
            parse_repeatmasker_out(p)

    def test_consensus_length_uses_mode_not_max(self, golden_dir, tmp_path):
        # one artifact hit implying an inflated model must not win over the mode
        src = open(golden_dir["out"]).read()
        inflated = (
            " 2001  1.0  0.0  0.0  chrT     95001    95200 (4800) + "
            " TST1       LINE/TST            201    400  (1600) 77\n"
        )  # implies model length 2000 vs the 7 consistent hits' 1000
        p = tmp_path / "inflated.out"
        p.write_text(src + inflated)
        df = parse_repeatmasker_out(p)
        assert consensus_lengths(df)["TST1"] == gd.CONS_LEN
