"""Local cache for database feature lookups (local-first, offline-capable).

Layout: ``<TESCOPE_CACHE_DIR>/features/dfam/<rm_name>__<cons_length>.json``

Each entry is the full lookup snapshot (features + provenance +
candidates) exactly as returned by the provider, so an old analysis can
be reproduced without any network access. The cache key includes the
expected consensus length because compatibility is evaluated against it.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from ..reference.store import default_cache_dir


def _slug(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name) or "unnamed"


def cache_path(rm_name: str, cons_length: int, provider: str = "dfam",
               cache_dir=None) -> Path:
    base = Path(cache_dir) if cache_dir else default_cache_dir()
    return base / "features" / provider / f"{_slug(rm_name)}__{int(cons_length)}.json"


def load_cached(rm_name: str, cons_length: int, provider: str = "dfam",
                cache_dir=None) -> dict | None:
    p = cache_path(rm_name, cons_length, provider, cache_dir)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None  # unreadable cache = no cache; provider may refresh it


def store_cached(rm_name: str, cons_length: int, snapshot: dict,
                 provider: str = "dfam", cache_dir=None) -> Path:
    p = cache_path(rm_name, cons_length, provider, cache_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(snapshot, indent=2), encoding="utf-8")
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
