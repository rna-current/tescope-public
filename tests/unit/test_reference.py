"""Reference store: custom import, caching, provenance, queries."""

from __future__ import annotations

import json

import golden_dataset as gd
import pytest

from tescope.core.annotation.frame_cache import load_frame, save_frame
from tescope.core.reference.registry import KNOWN_ASSEMBLIES
from tescope.core.reference.store import (
    ReferenceError,
    ReferenceNotInstalledError,
    ReferenceStore,
)


class TestFrameCache:
    def test_roundtrip(self, golden_dir, tmp_path):
        from tescope.core.annotation.repeatmasker import parse_repeatmasker_out

        df = parse_repeatmasker_out(golden_dir["out"])
        p = tmp_path / "f.npz"
        save_frame(df, p)
        back = load_frame(p)
        assert len(back) == len(df)
        assert list(back.columns) == list(df.columns)
        for col in df.columns:
            assert back[col].astype(str).tolist() == df[col].astype(str).tolist(), col


class TestStore:
    def test_custom_install_and_query(self, golden_store):
        assert golden_store.is_installed("testasm")
        assert "testasm" in golden_store.installed_assemblies()
        df = golden_store.load_annotation("testasm")
        assert len(df) == len(gd.COPIES) + 1
        assert golden_store.consensus_lengths("testasm")["TST1"] == gd.CONS_LEN
        sizes = golden_store.chrom_sizes("testasm")
        assert sizes == {gd.CHROM: gd.CHROM_LEN}

    def test_provenance_recorded(self, golden_store):
        prov = golden_store.provenance("testasm")
        assert prov.assembly == "testasm"
        assert prov.source == "user-file"
        assert prov.source_format == "repeatmasker-out"
        assert prov.md5
        assert prov.retrieved_at
        assert prov.extra["n_hits"] == len(gd.COPIES) + 1
        # original source file is retained
        src_dir = golden_store.assembly_dir("testasm") / "source"
        assert any(src_dir.iterdir())

    def test_select_repeat(self, golden_store):
        sub = golden_store.select_repeat("testasm", "TST1")
        assert len(sub) == len(gd.COPIES)

    def test_select_repeat_suggests(self, golden_store):
        with pytest.raises(ReferenceError, match="Did you mean.*TST1"):
            golden_store.select_repeat("testasm", "tst")

    def test_not_installed_message(self, tmp_path):
        store = ReferenceStore(cache_dir=tmp_path)
        with pytest.raises(ReferenceNotInstalledError, match="reference fetch hg38"):
            store.load_annotation("hg38")
        with pytest.raises(ReferenceNotInstalledError, match="reference add"):
            store.load_annotation("myasm")

    def test_double_install_requires_force(self, golden_store, golden_dir):
        with pytest.raises(ReferenceError, match="force"):
            golden_store.add_custom("testasm", golden_dir["out"])

    def test_hierarchy(self, golden_store):
        h = golden_store.hierarchy("testasm")
        assert set(h["repeat_name"].astype(str)) == {"TST1", "TST2"}

    def test_cache_rebuild_on_format_drift(self, golden_store):
        # corrupt the cache; load must transparently rebuild from source
        p = golden_store.assembly_dir("testasm") / "annotation.npz"
        p.write_bytes(b"garbage")
        golden_store._annotation_cache.clear()
        from tescope.core.annotation.frame_cache import FrameCacheError

        with pytest.raises(FrameCacheError):
            load_frame(p)
        df = golden_store.load_annotation("testasm")
        assert len(df) == len(gd.COPIES) + 1

    def test_known_assemblies_registry(self):
        assert {"hg38", "hg19", "mm39", "mm10"} <= set(KNOWN_ASSEMBLIES)
        for spec in KNOWN_ASSEMBLIES.values():
            assert spec.rmsk_url.startswith("https://")
            assert spec.chromsizes_url.startswith("https://")

    def test_fetch_unknown_assembly(self, tmp_path):
        store = ReferenceStore(cache_dir=tmp_path)
        with pytest.raises(ReferenceError, match="not a known assembly"):
            store.fetch("hg99")

    def test_provenance_json_is_valid(self, golden_store):
        p = golden_store.assembly_dir("testasm") / "provenance.json"
        data = json.loads(p.read_text())
        assert data["assembly"] == "testasm"
