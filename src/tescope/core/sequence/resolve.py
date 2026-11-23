"""Sequence resolution: config -> SequenceResolution, cache-first.

Kept strictly separate from projection/aggregation: sequence resolution
failure degrades to "unavailable" with an explanation and NEVER fails the
scientific analysis — sequence availability is never required to plot a
normal TEscope profile.
"""

from __future__ import annotations

from ..features.dfam import DfamClient, DfamError
from .cache import load_cached, store_cached
from .dfam import sequence_lookup
from .library import library_lookup
from .model import STATUS_UNAVAILABLE, SequenceResolution


def resolve_sequence(
    seq_cfg,
    repeat_name: str,
    cons_length: int,
    cache_dir=None,
    client: DfamClient | None = None,
    allow_network: bool = True,
) -> SequenceResolution:
    """Build the SequenceResolution for an analysis.

    ``seq_cfg`` is :class:`tescope.config.model.SequenceConfig`.
    Dfam is cache-first: a cached snapshot (verified or not) answers
    offline; the network is only consulted when nothing is cached and
    ``allow_network`` is true.
    """
    source = seq_cfg.source
    if source in ("", "none"):
        return SequenceResolution(
            status=STATUS_UNAVAILABLE, detail="no sequence source configured"
        )

    if source == "dfam":
        cached = load_cached(repeat_name, cons_length, "dfam", cache_dir)
        if cached is not None:
            cached.detail = f"{cached.detail} (cached)".strip()
            return cached
        if not allow_network:
            return SequenceResolution(
                status=STATUS_UNAVAILABLE,
                detail="Dfam sequence not cached and network lookups are disabled",
            )
        try:
            res = sequence_lookup(repeat_name, cons_length, client=client)
        except DfamError as e:
            return SequenceResolution(
                status=STATUS_UNAVAILABLE,
                detail=f"Dfam unreachable and no cached sequence ({e})",
            )
        store_cached(repeat_name, cons_length, res, "dfam", cache_dir)
        return res

    if source in ("library", "fasta"):
        if not seq_cfg.path:
            return SequenceResolution(
                status=STATUS_UNAVAILABLE,
                detail=f"sequence source '{source}' requires a file path",
            )
        provider = (
            "user-provided FASTA" if source == "fasta" else "user-provided library"
        )
        return library_lookup(
            seq_cfg.path,
            seq_cfg.family or repeat_name,
            cons_length,
            provider=provider,
            library_name=seq_cfg.library_name,
        )

    return SequenceResolution(
        status=STATUS_UNAVAILABLE,
        detail=f"unknown sequence source {source!r}",
    )
