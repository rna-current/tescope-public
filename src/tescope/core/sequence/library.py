"""User-provided consensus libraries (FASTA).

TEscope reads locally provided repeat libraries — a researcher's custom
consensus set, a RepeatMasker-associated library, or a licensed library
the user has independently and lawfully obtained. TEscope never bundles,
redistributes or automatically downloads licensed databases; it simply
reads the local file the user points it to (docs/sequence-browser.md,
"Licensing"). Provenance uses the neutral term "user-provided".

Supported: standard (optionally gzip-compressed) FASTA, including
RepeatMasker-style ``name#class/family`` headers. Family resolution is
careful and never guesses: an exact header-token match wins, then an
exact match on the part before ``#``; multiple matches are an error that
lists the candidates.

Large libraries are scanned once per (path, mtime, size) and the header
index is kept in memory for the process (no repeated re-parsing).
"""

from __future__ import annotations

import gzip
import os
from datetime import datetime, timezone
from pathlib import Path

from ..features.model import Compatibility, verify_length_compatibility
from .model import (
    STATUS_INCOMPATIBLE,
    STATUS_UNAVAILABLE,
    STATUS_VERIFIED,
    ConsensusSequence,
    SequenceError,
    SequenceProvenance,
    SequenceResolution,
    clean_sequence,
    sequence_checksum,
)

# process-lifetime index cache: (abspath, mtime_ns, size) -> records
_INDEX_CACHE: dict[tuple, dict[str, str]] = {}
_INDEX_CACHE_MAX = 8


def _open_text(path: Path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return open(path, encoding="utf-8", errors="strict")


def read_fasta_records(path: str | Path) -> dict[str, str]:
    """All records of a FASTA file: header first token -> raw sequence.

    Duplicate header tokens are an error (a library where one name maps to
    two different sequences cannot be resolved without guessing).
    """
    p = Path(path)
    if not p.is_file():
        raise SequenceError(f"library file not found: {p}")
    stat = os.stat(p)
    key = (str(p.resolve()), stat.st_mtime_ns, stat.st_size)
    cached = _INDEX_CACHE.get(key)
    if cached is not None:
        return cached

    records: dict[str, str] = {}
    name = None
    chunks: list[str] = []

    def flush():
        if name is None:
            return
        if name in records:
            raise SequenceError(
                f"{p.name}: duplicate FASTA record name {name!r} — "
                "cannot resolve families unambiguously"
            )
        records[name] = "".join(chunks)

    try:
        with _open_text(p) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    flush()
                    header = line[1:].strip()
                    if not header:
                        raise SequenceError(f"{p.name}: empty FASTA header")
                    name = header.split()[0]
                    chunks = []
                elif line.strip():
                    if name is None:
                        raise SequenceError(
                            f"{p.name}: sequence data before the first FASTA header"
                        )
                    chunks.append(line.strip())
        flush()
    except (OSError, UnicodeDecodeError) as e:
        raise SequenceError(f"cannot read library {p}: {e}") from e
    if not records:
        raise SequenceError(f"{p.name}: no FASTA records found")

    if len(_INDEX_CACHE) >= _INDEX_CACHE_MAX:
        _INDEX_CACHE.pop(next(iter(_INDEX_CACHE)))
    _INDEX_CACHE[key] = records
    return records


def _base_name(token: str) -> str:
    """RepeatMasker library headers are ``name#class/family``."""
    return token.split("#", 1)[0]


def find_family(records: dict[str, str], family: str) -> tuple[str, str]:
    """Resolve ``family`` to (header_token, matched_by) — never a guess.

    Exact token match first, then exact match on the pre-``#`` name.
    Zero matches raises with nearby names; multiple matches raises with
    all candidates listed.
    """
    if family in records:
        return family, "header"
    matches = [tok for tok in records if _base_name(tok) == family]
    if len(matches) == 1:
        return matches[0], "header name (RepeatMasker name#class/family)"
    if len(matches) > 1:
        raise SequenceError(
            f"family {family!r} is ambiguous in this library: "
            + ", ".join(sorted(matches)[:8])
        )
    near = sorted(
        tok for tok in records if family.lower() in _base_name(tok).lower()
    )[:5]
    hint = f" (similar: {', '.join(near)})" if near else ""
    raise SequenceError(f"family {family!r} not found in this library{hint}")


def library_lookup(
    path: str | Path,
    family: str,
    expected_cons_length: int,
    provider: str = "user-provided library",
    library_name: str = "",
) -> SequenceResolution:
    """Resolve one family from a user-provided FASTA library.

    A single-record file resolves without a family match when ``family``
    is empty or absent from the file's headers is NOT accepted silently —
    the single record is used only when explicitly the only one, and
    ``matched_by`` records that.
    """
    p = Path(path)
    base = SequenceProvenance(
        provider=provider,
        family=family,
        source_version=library_name,
        source_path=str(p),
        retrieved_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        expected_model_length=int(expected_cons_length),
    )
    try:
        records = read_fasta_records(p)
        if family:
            try:
                token, matched_by = find_family(records, family)
            except SequenceError:
                if len(records) == 1:
                    token, matched_by = next(iter(records)), "only record in file"
                else:
                    raise
        elif len(records) == 1:
            token, matched_by = next(iter(records)), "only record in file"
        else:
            raise SequenceError(
                f"{p.name} contains {len(records)} records; a family name is required"
            )
        seq = clean_sequence(records[token], f"{p.name}:{token}")
    except SequenceError as e:
        base.compatibility = STATUS_UNAVAILABLE
        base.compatibility_detail = str(e)
        return SequenceResolution(
            status=STATUS_UNAVAILABLE, detail=str(e), provenance=base
        )

    base.record_name = token
    base.matched_by = matched_by
    base.sequence_length = len(seq)
    base.checksum = sequence_checksum(seq)

    status, detail = verify_length_compatibility(expected_cons_length, len(seq))
    if status is Compatibility.VERIFIED:
        base.compatibility = STATUS_VERIFIED
        base.compatibility_detail = detail
        return SequenceResolution(
            status=STATUS_VERIFIED,
            detail=f"{token} ({p.name}) — {detail}",
            sequence=ConsensusSequence(sequence=seq, provenance=base),
            provenance=base,
        )
    base.compatibility = STATUS_INCOMPATIBLE
    base.compatibility_detail = detail
    return SequenceResolution(
        status=STATUS_INCOMPATIBLE,
        detail=f"{token} ({p.name}): {detail}",
        provenance=base,
    )
