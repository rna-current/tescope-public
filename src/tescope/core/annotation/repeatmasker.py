"""Parser for native RepeatMasker ``.out`` files.

Format reference (RepeatMasker documentation):

```
   SW  perc perc perc  query      position in query           matching  repeat         position in repeat
score  div. del. ins.  sequence    begin  end        (left)   repeat    class/family   begin  end  (left)  ID
```

- Two header lines followed by a blank line (some files omit headers).
- ``strand`` column is ``+`` or ``C`` (complement).
- Query coordinates are 1-based inclusive.
- Repeat (consensus) coordinates are 1-based inclusive and always given in
  the consensus orientation. For ``+`` hits the last three repeat fields
  are ``begin end (left)``; for ``C`` hits they are ``(left) end begin``.
- An optional trailing ``*`` marks a hit that is overlapped by a
  higher-scoring hit.

All coordinates are converted to 0-based half-open in the canonical frame.
"""

from __future__ import annotations

import gzip
import io
import os
from collections.abc import Iterable

import pandas as pd

from .schema import AnnotationSchemaError, finalize_frame

# class/family field is "class/family" or just "class"
_HEADER_TOKENS = ("SW", "score", "bit")


class RepeatMaskerParseError(AnnotationSchemaError):
    pass


def _open_text(path: str | os.PathLike) -> io.TextIOBase:
    p = os.fspath(path)
    if p.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")
    return open(p, encoding="utf-8", errors="replace")


def _strip_parens(tok: str) -> tuple[int, bool]:
    """Return (value, was_parenthesized)."""
    if tok.startswith("(") and tok.endswith(")"):
        return int(tok[1:-1]), True
    return int(tok), False


def _split_class_family(tok: str) -> tuple[str, str]:
    if "/" in tok:
        cls, fam = tok.split("/", 1)
        return cls, fam
    return tok, tok


def parse_repeatmasker_out(
    path: str | os.PathLike,
    *,
    keep_secondary: bool = True,
) -> pd.DataFrame:
    """Parse a RepeatMasker ``.out`` (or ``.out.gz``) file to the canonical frame.

    Parameters
    ----------
    path:
        Path to the ``.out`` file, optionally gzip-compressed.
    keep_secondary:
        Keep hits flagged with a trailing ``*`` (overlapped by a
        higher-scoring hit). They are retained by default and marked in the
        ``is_secondary`` column so filtering stays explicit.
    """
    rows = []
    n_bad = 0
    first_bad: tuple[int, str] | None = None
    with _open_text(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            fields = line.split()
            # Header lines start with "SW"/"score"/"perc" style tokens.
            if fields[0] in _HEADER_TOKENS or fields[0] in ("perc", "There"):
                continue
            try:
                row = _parse_line(fields)
            except (ValueError, IndexError):
                # Real RepeatMasker output contains rare artifact lines (e.g.
                # consensus begin == end + 1 on tiny hits). Skip and count;
                # only a large malformed fraction indicates the wrong format.
                n_bad += 1
                if first_bad is None:
                    first_bad = (lineno, line[:200])
                if n_bad > 10_000 and n_bad > len(rows):
                    raise RepeatMaskerParseError(
                        f"{os.fspath(path)}: mostly unparseable as RepeatMasker .out; "
                        f"first bad line {first_bad[0]}: {first_bad[1]!r}"
                    ) from None
                continue
            rows.append(row)
    if not rows:
        raise RepeatMaskerParseError(
            f"{os.fspath(path)}: no RepeatMasker hits could be parsed"
            + (f" ({n_bad} malformed lines)" if n_bad else "")
        )
    if n_bad > 100 and n_bad > 0.05 * (n_bad + len(rows)):
        raise RepeatMaskerParseError(
            f"{os.fspath(path)}: {n_bad} of {n_bad + len(rows)} data lines "
            f"unparseable — not valid RepeatMasker .out? "
            f"first bad line {first_bad[0]}: {first_bad[1]!r}"
        )
    df = pd.DataFrame(
        rows,
        columns=[
            "sw_score",
            "pct_div",
            "pct_del",
            "pct_ins",
            "chrom",
            "genome_start",
            "genome_end",
            "strand",
            "repeat_name",
            "repeat_class",
            "repeat_family",
            "cons_start",
            "cons_end",
            "cons_left",
            "rm_id",
            "is_secondary",
        ],
    )
    if not keep_secondary:
        df = df[~df["is_secondary"]]
    out = finalize_frame(df)
    out.attrs["n_malformed_lines"] = n_bad
    return out


def _parse_line(f: list[str]) -> tuple:
    """Parse one whitespace-split .out data line."""
    sw_score = float(f[0])
    pct_div = float(f[1])
    pct_del = float(f[2])
    pct_ins = float(f[3])
    chrom = f[4]
    q_begin = int(f[5])
    q_end = int(f[6])
    # f[7] is "(left)" for the query; unused.
    strand_tok = f[8]
    if strand_tok not in ("+", "C"):
        raise ValueError(f"bad strand token {strand_tok!r}")
    strand = "+" if strand_tok == "+" else "-"
    repeat_name = f[9]
    repeat_class, repeat_family = _split_class_family(f[10])

    r1, p1 = _strip_parens(f[11])
    r2, _ = _strip_parens(f[12])
    r3, p3 = _strip_parens(f[13])
    if strand == "+":
        # begin end (left)
        if p1 or not p3:
            raise ValueError("unexpected parenthesization for + strand repeat coords")
        cons_begin, cons_end_incl, cons_left = r1, r2, r3
    else:
        # (left) end begin
        if not p1 or p3:
            raise ValueError("unexpected parenthesization for C strand repeat coords")
        cons_left, cons_end_incl, cons_begin = r1, r2, r3
    if cons_begin < 1 or cons_end_incl < cons_begin:
        raise ValueError(f"bad consensus coords {cons_begin}..{cons_end_incl}")

    rm_id = -1
    is_secondary = False
    rest = f[14:]
    if rest:
        if rest[-1] == "*":
            is_secondary = True
            rest = rest[:-1]
        if rest:
            try:
                rm_id = int(rest[0])
            except ValueError:
                rm_id = -1

    return (
        sw_score,
        pct_div,
        pct_del,
        pct_ins,
        chrom,
        q_begin - 1,  # 1-based inclusive -> 0-based half-open
        q_end,
        strand,
        repeat_name,
        repeat_class,
        repeat_family,
        cons_begin - 1,
        cons_end_incl,
        cons_left,
        rm_id,
        is_secondary,
    )


def sniff_format(path: str | os.PathLike) -> str:
    """Best-effort format detection: 'repeatmasker-out' or 'ucsc-rmsk'.

    UCSC rmsk dumps are 17-column tab-separated with a numeric ``bin``
    first column and strand at column 10; .out files are whitespace-aligned
    with strand ('+'/'C') at column 9.
    """
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) >= 16 and fields[9] in ("+", "-"):
                return "ucsc-rmsk"
            wfields = line.split()
            if wfields and wfields[0] in _HEADER_TOKENS or (wfields and wfields[0] == "perc"):
                return "repeatmasker-out"
            if len(wfields) >= 14 and wfields[8] in ("+", "C"):
                return "repeatmasker-out"
            # keep scanning a few lines
    raise RepeatMaskerParseError(
        f"{os.fspath(path)}: could not detect RepeatMasker format "
        "(expected native .out or UCSC rmsk table)"
    )


def parse_any(path: str | os.PathLike) -> tuple[pd.DataFrame, str]:
    """Parse either supported format; returns (frame, detected_format)."""
    fmt = sniff_format(path)
    if fmt == "ucsc-rmsk":
        from .ucsc_rmsk import parse_ucsc_rmsk

        return parse_ucsc_rmsk(path), fmt
    return parse_repeatmasker_out(path), fmt


def iter_paths(paths: Iterable[str | os.PathLike]) -> pd.DataFrame:
    """Parse and concatenate several annotation files (same assembly)."""
    frames = [parse_any(p)[0] for p in paths]
    if len(frames) == 1:
        return frames[0]
    return finalize_frame(pd.concat(frames, ignore_index=True))
