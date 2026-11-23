from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))  # make golden_dataset importable

import golden_dataset  # noqa: E402


@pytest.fixture(scope="session")
def golden_dir(tmp_path_factory) -> dict[str, str]:
    """Paths of the materialized golden dataset (built once per session)."""
    d = tmp_path_factory.mktemp("golden")
    return golden_dataset.build_all(d)


@pytest.fixture(scope="session")
def golden_store(golden_dir, tmp_path_factory):
    """A ReferenceStore with the golden annotation installed as 'testasm'."""
    from tescope.core.reference.store import ReferenceStore

    cache = tmp_path_factory.mktemp("cache")
    store = ReferenceStore(cache_dir=cache)
    store.add_custom(
        "testasm",
        golden_dir["out"],
        chrom_sizes_path=golden_dir["chrom_sizes"],
        notes="synthetic golden dataset",
    )
    return store
