"""Dfam consensus-sequence provider.

Reuses the existing Dfam integration (:mod:`tescope.core.features.dfam`):
the same official API, the same name/alias matching, and the same
exact-length compatibility gate. Only a VERIFIED candidate's sequence is
fetched; when name-matching candidates exist but none is length-
compatible (e.g. current Dfam represents L1HS as split 5'/3' fragments
while the RepeatMasker model is full-length), the result is INCOMPATIBLE
with every rejection recorded — TEscope never concatenates fragments or
substitutes a different model to force a fit.

Dfam data is released under CC0 1.0, so retrieving and caching sequence
records is permitted.
"""

from __future__ import annotations

from datetime import datetime, timezone

from ..features.dfam import DfamClient, DfamError, _matches_rm_name
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


def _parse_single_fasta(text: str, where: str) -> str:
    """The sequence of a single-record FASTA response."""
    lines = text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise SequenceError(f"{where}: response is not FASTA")
    body = []
    for line in lines[1:]:
        if line.startswith(">"):
            break  # single-family endpoint: only the first record is ours
        body.append(line.strip())
    return clean_sequence("".join(body), where)


def sequence_lookup(
    rm_name: str,
    expected_cons_length: int,
    client: DfamClient | None = None,
) -> SequenceResolution:
    """Look up the Dfam consensus sequence for a RepeatMasker repeat name.

    Same candidate discovery as the feature lookup (name match, then
    alias check on full records — UCSC RepeatMasker uses Repbase names,
    which Dfam carries as aliases). Raises :class:`DfamError` on network
    failure so callers can fall back to cache/degrade gracefully.
    """
    client = client or DfamClient()
    version = client.version()
    db_version = f"Dfam {version.get('dfam_version', '?')}"
    release_date = str(version.get("release_date", ""))
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")

    results = client.search_families(rm_name)
    candidates = [r for r in results if _matches_rm_name(r, rm_name)]
    if not candidates:
        for r in results[:10]:
            try:
                full = client.family(r["accession"])
            except DfamError:
                continue
            if _matches_rm_name(full, rm_name):
                full.setdefault("length", r.get("length"))
                candidates.append(full)

    base = SequenceProvenance(
        provider="Dfam",
        family=rm_name,
        source_version=db_version,
        source_release_date=release_date,
        retrieved_at=now,
        api_url=client.base_url,
        expected_model_length=int(expected_cons_length),
    )

    if not candidates:
        base.compatibility = STATUS_UNAVAILABLE
        base.compatibility_detail = (
            f"no Dfam family with name or alias '{rm_name}' ({db_version})"
        )
        return SequenceResolution(
            status=STATUS_UNAVAILABLE,
            detail=base.compatibility_detail,
            provenance=base,
        )

    rejections = []
    for r in candidates:
        status, detail = verify_length_compatibility(
            expected_cons_length, r.get("length")
        )
        if status is not Compatibility.VERIFIED:
            rejections.append(f"{r.get('accession')} ({r.get('name')}): {detail}")
            continue
        seq = _parse_single_fasta(
            client.family_sequence_fasta(r["accession"]),
            f"Dfam sequence {r.get('accession')}",
        )
        prov = base
        prov.accession = str(r.get("accession", ""))
        prov.record_name = str(r.get("name", ""))
        prov.matched_by = _matches_rm_name(r, rm_name) or "name"
        prov.sequence_length = len(seq)
        prov.checksum = sequence_checksum(seq)
        # honesty check: the record SAID it matched; the actual sequence
        # must too, or the coordinate promise would be false
        if len(seq) != int(expected_cons_length):
            prov.compatibility = STATUS_INCOMPATIBLE
            prov.compatibility_detail = (
                f"{prov.accession}: record length matched but the fetched "
                f"sequence is {len(seq)} bp (expected "
                f"{expected_cons_length} bp)"
            )
            rejections.append(prov.compatibility_detail)
            continue
        prov.compatibility = STATUS_VERIFIED
        prov.compatibility_detail = detail
        return SequenceResolution(
            status=STATUS_VERIFIED,
            detail=f"{prov.accession} '{prov.record_name}' — {detail}",
            sequence=ConsensusSequence(sequence=seq, provenance=prov),
            provenance=prov,
        )

    base.compatibility = STATUS_INCOMPATIBLE
    base.compatibility_detail = "; ".join(rejections)
    return SequenceResolution(
        status=STATUS_INCOMPATIBLE,
        detail=base.compatibility_detail,
        provenance=base,
    )
