"""Consensus-feature annotation subsystem tests.

All Dfam interactions are mocked with clearly SYNTHETIC fixture records
(names like SYNTH_*): no biological feature name or coordinate here is a
real-world claim, and no network access occurs. A live-API check exists
separately, gated behind TESCOPE_LIVE_DFAM=1.
"""

from __future__ import annotations

import json
import os

import pytest

from tescope.config.model import AnalysisConfig, AnnotationConfig, RepeatSelection, TrackConfig
from tescope.core.features.cache import cache_path, drop_cached, load_cached, store_cached
from tescope.core.features.dfam import DfamClient, DfamError, lookup
from tescope.core.features.model import (
    SOURCE_DATABASE,
    SOURCE_USER,
    Compatibility,
    ConsensusFeature,
    FeatureSet,
    verify_length_compatibility,
)
from tescope.core.features.resolve import AnnotationRequiredError, resolve_features
from tescope.core.features.user import UserFeatureError, parse_user_features_file


class FakeDfamClient(DfamClient):
    """Offline stand-in serving synthetic records."""

    def __init__(self, families: list[dict], fail: bool = False):
        super().__init__(base_url="https://synthetic.invalid/api")
        self._families = {f["accession"]: f for f in families}
        self._fail = fail
        self.calls: list[str] = []

    def _get_json(self, path, params=None):
        self.calls.append(path)
        if self._fail:
            raise DfamError("synthetic network failure")
        if path == "/version":
            return {"dfam_version": "9.9-synthetic", "release_date": "TEST"}
        if path == "/families":
            name = params["name"]
            return {"results": [
                {k: f.get(k) for k in ("accession", "name", "length")}
                for f in self._families.values()
                if name in f.get("name", "")
                or any(a["alias"] == name for a in f.get("aliases", []))
            ]}
        acc = path.rsplit("/", 1)[-1]
        if acc in self._families:
            return self._families[acc]
        raise DfamError(f"404 {path}")


SYNTH_FAMILY = {
    "accession": "DFTEST000001",
    "name": "SYNTH_FAM",
    "version": 7,
    "length": 1000,
    "aliases": [{"database": "Repbase", "alias": "TST1"}],
    "features": [
        {
            "type": "protein_match",
            "description": "synthetic evidence record",
            "model_start_pos": 101,       # 1-based inclusive (Dfam convention)
            "model_end_pos": 400,
            "label": "SYNTH_MATCH_A",
            "attributes": [{"attribute": "evalue", "value": "1e-10"}],
        },
    ],
    "coding_seqs": [
        {"product": "SYNTH_CDS_B", "cds_start": 500, "cds_end": 900,
         "description": "synthetic cds"},
    ],
}


def make_client(length=1000, fail=False):
    fam = json.loads(json.dumps(SYNTH_FAMILY))
    fam["length"] = length
    return FakeDfamClient([fam], fail=fail)


class TestCompatibility:
    def test_exact_match_verified(self):
        st, detail = verify_length_compatibility(1000, 1000)
        assert st is Compatibility.VERIFIED and "1000" in detail

    def test_mismatch_incompatible(self):
        st, detail = verify_length_compatibility(1000, 902)
        assert st is Compatibility.INCOMPATIBLE and "mismatch" in detail

    def test_missing_unverified(self):
        st, _ = verify_length_compatibility(1000, None)
        assert st is Compatibility.UNVERIFIED


class TestDfamLookup:
    def test_verified_lookup_extracts_verbatim_features(self):
        snap = lookup("TST1", 1000, client=make_client())
        prov = snap["provenance"]
        assert prov["compatibility"] == "verified"
        assert prov["accession"] == "DFTEST000001"
        assert prov["database_version"] == "Dfam 9.9-synthetic"
        assert prov["matched_by"].startswith("alias:")
        names = [f["name"] for f in snap["features"]]
        assert names == ["SYNTH_MATCH_A", "SYNTH_CDS_B"]
        # 1-based inclusive -> 0-based half-open; types passed through
        fa = snap["features"][0]
        assert (fa["start"], fa["end"]) == (100, 400)
        assert fa["type"] == "protein_match"  # evidence terminology preserved
        assert fa["source_attributes"]["evalue"] == "1e-10"
        cds = snap["features"][1]
        assert (cds["start"], cds["end"]) == (499, 900)

    def test_length_mismatch_withholds_features(self):
        snap = lookup("TST1", 1200, client=make_client(length=1000))
        assert snap["provenance"]["compatibility"] == "incompatible"
        assert snap["features"] == []
        assert "mismatch" in snap["provenance"]["compatibility_detail"]

    def test_unknown_name_unverified(self):
        snap = lookup("NOSUCH", 1000, client=make_client())
        assert snap["provenance"]["compatibility"] == "unverified"
        assert snap["features"] == []


class TestCache:
    def test_roundtrip_and_refresh(self, tmp_path):
        snap = lookup("TST1", 1000, client=make_client())
        p = store_cached("TST1", 1000, snap, cache_dir=tmp_path)
        assert p == cache_path("TST1", 1000, cache_dir=tmp_path)
        back = load_cached("TST1", 1000, cache_dir=tmp_path)
        assert back == snap
        assert drop_cached("TST1", 1000, cache_dir=tmp_path)
        assert load_cached("TST1", 1000, cache_dir=tmp_path) is None

    def test_key_includes_consensus_length(self, tmp_path):
        store_cached("TST1", 1000, {"x": 1}, cache_dir=tmp_path)
        assert load_cached("TST1", 999, cache_dir=tmp_path) is None


class TestResolve:
    def cfg(self, **kw):
        return AnnotationConfig(enabled=True, provider="dfam", **kw)

    def test_network_then_cache(self, tmp_path):
        client = make_client()
        fs = resolve_features(self.cfg(), "TST1", 1000, cache_dir=tmp_path,
                              client=client)
        assert len(fs.database_features) == 2
        assert fs.provenance.compatibility == "verified"
        calls_after_first = len(client.calls)
        # second resolve: cache only, no client calls
        fs2 = resolve_features(self.cfg(), "TST1", 1000, cache_dir=tmp_path,
                               client=client)
        assert len(client.calls) == calls_after_first
        assert fs2.database_features[0].name == "SYNTH_MATCH_A"
        assert "cached" in fs2.status

    def test_offline_with_cache(self, tmp_path):
        resolve_features(self.cfg(), "TST1", 1000, cache_dir=tmp_path,
                         client=make_client())
        fs = resolve_features(self.cfg(), "TST1", 1000, cache_dir=tmp_path,
                              client=make_client(fail=True), allow_network=False)
        assert len(fs.database_features) == 2  # served from cache, no network

    def test_network_failure_degrades_gracefully(self, tmp_path):
        fs = resolve_features(self.cfg(), "TST1", 1000, cache_dir=tmp_path,
                              client=make_client(fail=True))
        assert fs.database_features == []
        assert "unavailable" in fs.status

    def test_required_raises_only_when_nothing(self, tmp_path):
        with pytest.raises(AnnotationRequiredError):
            resolve_features(self.cfg(required=True), "TST1", 1000,
                             cache_dir=tmp_path, client=make_client(fail=True))
        fs = resolve_features(
            self.cfg(required=True,
                     user_features=[{"name": "u", "start": 0, "end": 10}]),
            "TST1", 1000, cache_dir=tmp_path, client=make_client(fail=True))
        assert len(fs.user_features) == 1

    def test_incompatible_never_displayed(self, tmp_path):
        fs = resolve_features(self.cfg(), "TST1", 1200, cache_dir=tmp_path,
                              client=make_client(length=1000))
        assert fs.database_features == []
        assert "incompatible" in fs.status

    def test_origins_never_conflated(self, tmp_path):
        fs = resolve_features(
            self.cfg(user_features=[{"name": "mine", "start": 5, "end": 50}]),
            "TST1", 1000, cache_dir=tmp_path, client=make_client())
        assert all(f.origin == SOURCE_DATABASE for f in fs.database_features)
        assert all(f.origin == SOURCE_USER for f in fs.user_features)
        with pytest.raises(ValueError):
            bad = FeatureSet(repeat_name="x", cons_length=10,
                             database_features=[ConsensusFeature(
                                 name="u", start=0, end=5, origin=SOURCE_USER)])
            bad.validate()


class TestUserFeatures:
    def test_tsv_parse(self, tmp_path):
        p = tmp_path / "f.tsv"
        p.write_text("name\tstart\tend\tlabel\nregion1\t10\t100\tdemo\nregion2\t50\t70\n")
        feats = parse_user_features_file(p)
        assert [f.name for f in feats] == ["region1", "region2"]
        assert feats[0].origin == SOURCE_USER
        assert feats[0].type == "demo"

    def test_yaml_parse(self, tmp_path):
        p = tmp_path / "f.yaml"
        p.write_text("- name: a\n  start: 0\n  end: 10\n- {name: b, start: 5, end: 8, note: n}\n")
        feats = parse_user_features_file(p)
        assert len(feats) == 2 and feats[1].note == "n"

    @pytest.mark.parametrize("content,msg", [
        ("name\tstart\tend\nx\t10\t5\n", "must be >"),
        ("name\tstart\tend\nx\tten\t20\n", "integers"),
        ("start\tend\n1\t2\n", "header"),
        ("name\tstart\tend\n\t1\t2\n", "missing required"),
    ])
    def test_malformed_tsv(self, tmp_path, content, msg):
        p = tmp_path / "bad.tsv"
        p.write_text(content)
        with pytest.raises(UserFeatureError, match=msg):
            parse_user_features_file(p)

    def test_bounds_checked_at_attach(self, tmp_path):
        cfg = AnnotationConfig(enabled=True, provider="none",
                               user_features=[{"name": "big", "start": 0, "end": 5000}])
        with pytest.raises(ValueError, match="exceeds consensus"):
            resolve_features(cfg, "TST1", 1000, cache_dir=tmp_path)


class TestConfigRoundtrip:
    def test_yaml_roundtrip(self, tmp_path):
        from tescope.config.io import load_config, save_config

        cfg = AnalysisConfig(
            assembly="testasm", repeat=RepeatSelection(name="TST1"),
            tracks=[TrackConfig(file=str(tmp_path / "x.bw"))],
            annotations=AnnotationConfig(
                enabled=True, provider="none",
                user_features=[{"name": "u1", "start": 3, "end": 9, "label": "L"}],
                show_in_plots=False),
        )
        p = tmp_path / "a.yaml"
        save_config(cfg, p)
        cfg2 = load_config(p)
        assert cfg2.annotations.enabled is True
        assert cfg2.annotations.show_in_plots is False
        assert cfg2.annotations.user_features[0]["name"] == "u1"


@pytest.mark.skipif(os.environ.get("TESCOPE_LIVE_DFAM") != "1",
                    reason="live Dfam API check (set TESCOPE_LIVE_DFAM=1)")
class TestLiveDfam:
    def test_version_endpoint(self):
        v = DfamClient().version()
        assert "dfam_version" in v
