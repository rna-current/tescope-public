"""Operations on VERIFIED consensus sequences: subrange, FASTA, search.

All functions take a :class:`ConsensusSequence` that resolution has
already verified against the analysis model — callers must gate on
``SequenceResolution.verified`` (the GUI/CLI capability layer does).
Coordinates in and out are 1-based inclusive (see coords.py).

Search is LITERAL: IUPAC letters only, case-insensitive, no motif
prediction, no scoring, no biology — coordinate hits only. ``N`` in the
query matches any base (documented wildcard); other letters match
themselves exactly.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .coords import RangeError, base_range_to_slice
from .model import (
    ConsensusSequence,
    SequenceError,
    clean_sequence,
    reverse_complement,
)

_FASTA_WIDTH = 70


def subsequence(seq: ConsensusSequence, base_start: int, base_end: int) -> str:
    """Bases base_start..base_end (1-based inclusive)."""
    if base_start < 1 or base_end > seq.length or base_end < base_start:
        raise RangeError(
            f"range {base_start}-{base_end} outside 1-{seq.length} (1-based)"
        )
    return seq.sequence[base_range_to_slice(base_start, base_end)]


def fasta_text(
    seq: ConsensusSequence,
    base_start: int | None = None,
    base_end: int | None = None,
    tescope_version: str = "",
) -> str:
    """FASTA export with provenance in the header (explicit user action).

    Header records family, provider, accession/release where available,
    the 1-based coordinate range for partial exports, checksum and the
    TEscope version. User-provided sources keep their local path label;
    no redistribution rights are claimed.
    """
    b1 = 1 if base_start is None else base_start
    b2 = seq.length if base_end is None else base_end
    body = subsequence(seq, b1, b2)
    pv = seq.provenance
    parts = [pv.family or pv.record_name or "consensus"]
    src = pv.provider
    if pv.accession:
        src += f":{pv.accession}"
    if src:
        parts.append(f"source={src.replace(' ', '_')}")
    if pv.source_version:
        parts.append(f"release={pv.source_version.replace(' ', '_')}")
    if pv.source_path:
        parts.append(f"file={pv.source_path.replace(' ', '_')}")
    parts.append(f"positions={b1}-{b2}(1-based)")
    parts.append(f"length={b2 - b1 + 1}")
    if pv.checksum:
        parts.append(f"source_sha256={pv.checksum[:16]}")
    if tescope_version:
        parts.append(f"tescope={tescope_version}")
    lines = [">" + " ".join(parts)]
    for i in range(0, len(body), _FASTA_WIDTH):
        lines.append(body[i : i + _FASTA_WIDTH])
    return "\n".join(lines) + "\n"


@dataclass(frozen=True)
class SequenceHit:
    """One literal match: 1-based inclusive coordinates on the consensus."""

    base_start: int
    base_end: int
    strand: str  # "+" forward, "-" reverse complement of the query

    @property
    def label(self) -> str:
        s = f"{self.base_start}-{self.base_end}"
        return s if self.strand == "+" else f"{s} (revcomp)"


def search_literal(
    seq: ConsensusSequence,
    query: str,
    include_revcomp: bool = False,
    max_hits: int = 500,
) -> list[SequenceHit]:
    """All literal occurrences of ``query`` (1-based inclusive hits).

    Case-insensitive; ``N`` in the query matches any base. With
    ``include_revcomp``, occurrences of the reverse complement of the
    query are also reported, clearly flagged. Overlapping matches count.
    """
    q = clean_sequence(query, "search query")
    if not q:
        raise SequenceError("search query is empty")
    if len(q) > seq.length:
        return []

    def find(pattern: str, strand: str) -> list[SequenceHit]:
        rx = re.compile("(?=(" + re.escape(pattern).replace("N", ".") + "))")
        out = []
        for m in rx.finditer(seq.sequence):
            out.append(
                SequenceHit(m.start() + 1, m.start() + len(pattern), strand)
            )
            if len(out) >= max_hits:
                break
        return out

    hits = find(q, "+")
    if include_revcomp:
        rc = reverse_complement(q)
        if rc != q:
            hits += find(rc, "-")
    hits.sort(key=lambda h: (h.base_start, h.strand))
    return hits[:max_hits]
