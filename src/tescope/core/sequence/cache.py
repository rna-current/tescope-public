"""Local cache for Dfam sequence lookups (local-first, offline-capable).

Layout: ``<TESCOPE_CACHE_DIR>/sequences/dfam/<rm_name>__<cons_length>.json``

Each entry is a full :class:`SequenceResolution` snapshot (sequence +
provenance, including checksum and source release), so an old analysis
reproduces without network access. The cache key includes the expected
consensus model length — the identity is (family, coordinate model),
never the family name alone, so a model change in the annotation can
never silently reuse a sequence cached for a different model. The stored
provenance additionally records accession/release/checksum for auditing.

User-provided libraries are NOT copied into the cache: they are read
from their original location (the file stays under the user's control
and license).
"""

from __future__ import annotations

import json
from pathlib import Path

from ..features.cache import _slug
from ..reference.store import default_cache_dir
from .model import SequenceResolution


def cache_path(rm_name: str, cons_length: int, provider: str = "dfam",
               cache_dir=None) -> Path:
    base = Path(cache_dir) if cache_dir else default_cache_dir()
    return base / "sequences" / provider / f"{_slug(rm_name)}__{int(cons_length)}.json"


def load_cached(rm_name: str, cons_length: int, provider: str = "dfam",
                cache_dir=None) -> SequenceResolution | None:
    p = cache_path(rm_name, cons_length, provider, cache_dir)
    if not p.exists():
        return None
    try:
        return SequenceResolution.from_dict(
            json.loads(p.read_text(encoding="utf-8"))
        )
    except (OSError, ValueError, KeyError):
        return None  # unreadable cache = no cache; provider may refresh it


def store_cached(rm_name: str, cons_length: int, res: SequenceResolution,
                 provider: str = "dfam", cache_dir=None) -> Path:
    p = cache_path(rm_name, cons_length, provider, cache_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(res.to_dict(), indent=2), encoding="utf-8")
    tmp.replace(p)
    return p


def drop_cached(rm_name: str, cons_length: int, provider: str = "dfam",
                cache_dir=None) -> bool:
    """Remove one cache entry (refresh mechanism). Returns True if removed."""
    p = cache_path(rm_name, cons_length, provider, cache_dir)
    if p.exists():
        p.unlink()
        return True
    return False
