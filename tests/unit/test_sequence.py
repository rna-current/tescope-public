"""Consensus-sequence layer: providers, cache, compatibility, coordinates.

All sequences here are synthetic (no proprietary biological sequences in
fixtures). Coordinate tests pin the 0-based internal / 1-based displayed
convention with no off-by-one ambiguity.
"""

from __future__ import annotations

import gzip
import json

import pytest

from tescope.config.model import AnalysisConfig, ConfigError, SequenceConfig
from tescope.core.sequence.cache import (
    cache_path,
    drop_cached,
    load_cached,
    store_cached,
)
from tescope.core.sequence.capabilities import sequence_capabilities
from tescope.core.sequence.coords import (
    RangeError,
    axis_to_base,
    base_range_to_slice,
    base_to_axis_span,
    parse_position_or_range,
    visible_base_range,
)
from tescope.core.sequence.dfam import sequence_lookup
from tescope.core.sequence.library import (
    find_family,
    library_lookup,
    read_fasta_records,
)
from tescope.core.sequence.model import (
    STATUS_INCOMPATIBLE,
    STATUS_UNAVAILABLE,
    STATUS_VERIFIED,
    ConsensusSequence,
    SequenceError,
    SequenceProvenance,
    SequenceResolution,
    clean_sequence,
    reverse_complement,
    sequence_checksum,
)
from tescope.core.sequence.ops import fasta_text, search_literal, subsequence
from tescope.core.sequence.resolve import resolve_sequence

# a synthetic 40-base consensus with a distinct pattern (position -> base
# is unambiguous: base b (1-based) is SEQ40[b-1])
SEQ40 = "ACGTACGTGGGGCCCCAAAATTTTNNNNACGTTGCATGCA"


# ------------------------------------------------------------------ model
class TestModel:
    def test_clean_sequence_valid_iupac(self):
        assert clean_sequence("acgt\nRY SWkm bdhvn u") == "ACGTRYSWKMBDHVNU"

    def test_clean_sequence_invalid_character_position(self):
        with pytest.raises(SequenceError, match="base 5"):
            clean_sequence("ACGTX")

    def test_checksum_deterministic_and_distinct(self):
        assert sequence_checksum("ACGT") == sequence_checksum("ACGT")
        assert sequence_checksum("ACGT") != sequence_checksum("ACGA")

    def test_reverse_complement(self):
        assert reverse_complement("ACGTN") == "NACGT"
        assert reverse_complement(reverse_complement(SEQ40)) == SEQ40

    def test_resolution_roundtrip(self):
        res = SequenceResolution(
            status=STATUS_VERIFIED,
            detail="d",
            sequence=ConsensusSequence(
                SEQ40, SequenceProvenance(provider="p", checksum="c")),
            provenance=SequenceProvenance(provider="p", checksum="c"),
        )
        back = SequenceResolution.from_dict(
            json.loads(json.dumps(res.to_dict())))
        assert back.verified
        assert back.sequence.sequence == SEQ40
        assert back.provenance.provider == "p"


# ------------------------------------------------------------- coordinates
class TestCoordinates:
    """Internal 0-based half-open <-> displayed 1-based inclusive."""

    def test_axis_to_base_boundaries(self):
        assert axis_to_base(0.0, 40) == 1  # first base
        assert axis_to_base(0.999, 40) == 1
        assert axis_to_base(39.0, 40) == 40  # last base
        assert axis_to_base(39.999, 40) == 40
        assert axis_to_base(40.0, 40) is None  # past the end: not a base
        assert axis_to_base(-0.1, 40) is None

    def test_base_to_axis_span(self):
        # base 1 occupies [0, 1); bases 1..40 occupy [0, 40)
        assert base_to_axis_span(1, 1) == (0.0, 1.0)
        assert base_to_axis_span(1, 40) == (0.0, 40.0)
        assert base_to_axis_span(4000, 4500) == (3999.0, 4500.0)

    def test_slice_matches_display(self):
        seq = ConsensusSequence(SEQ40)
        for b in (1, 2, 39, 40):
            assert SEQ40[base_range_to_slice(b, b)] == SEQ40[b - 1]
            assert subsequence(seq, b, b) == SEQ40[b - 1]

    def test_visible_base_range(self):
        assert visible_base_range(0.0, 40.0, 40) == (1, 40)
        assert visible_base_range(3.2, 7.8, 40) == (4, 8)
        assert visible_base_range(-5.0, 100.0, 40) == (1, 40)

    def test_parse_position_and_range(self):
        assert parse_position_or_range("17", 40) == (17, 17)
        assert parse_position_or_range("4-8", 40) == (4, 8)
        assert parse_position_or_range(" 8 - 4 ", 40) == (4, 8)
        assert parse_position_or_range("1-40", 40) == (1, 40)

    def test_parse_rejects_out_of_bounds_and_garbage(self):
        for bad in ("0", "41", "1-41", "", "abc", "3--5"):
            with pytest.raises(RangeError):
                parse_position_or_range(bad, 40)


# ---------------------------------------------------------------- user FASTA
class TestUserLibrary:
    def _write(self, tmp_path, text, name="lib.fa"):
        p = tmp_path / name
        p.write_text(text)
        return p

    def test_multi_record_and_repeatmasker_headers(self, tmp_path):
        p = self._write(
            tmp_path,
            ">L1demo#LINE/L1 a demo\nACGTAC\nGTAC\n>Alu9#SINE/Alu\nGGGG\n",
        )
        recs = read_fasta_records(p)
        assert set(recs) == {"L1demo#LINE/L1", "Alu9#SINE/Alu"}
        assert recs["L1demo#LINE/L1"] == "ACGTACGTAC"
        tok, how = find_family(recs, "L1demo")
        assert tok == "L1demo#LINE/L1" and "name#class" in how

    def test_exact_header_beats_base_name(self, tmp_path):
        p = self._write(tmp_path, ">TST1\nAAAA\n>TST1#LINE/L1\nCCCC\n")
        recs = read_fasta_records(p)
        assert find_family(recs, "TST1")[0] == "TST1"

    def test_ambiguous_family_never_guessed(self, tmp_path):
        p = self._write(tmp_path, ">X#LINE/L1\nAAAA\n>X#SINE/Alu\nCCCC\n")
        with pytest.raises(SequenceError, match="ambiguous"):
            find_family(read_fasta_records(p), "X")

    def test_missing_family_lists_similar(self, tmp_path):
        p = self._write(tmp_path, ">L1demo#LINE/L1\nAAAA\n")
        with pytest.raises(SequenceError, match="similar: L1demo"):
            find_family(read_fasta_records(p), "l1de")

    def test_gzip_supported(self, tmp_path):
        p = tmp_path / "lib.fa.gz"
        with gzip.open(p, "wt") as fh:
            fh.write(">A\nACGT\n")
        assert read_fasta_records(p) == {"A": "ACGT"}

    def test_invalid_fasta_shapes(self, tmp_path):
        for text, msg in (
            ("ACGT\n", "before the first"),
            (">\nACGT\n", "empty FASTA header"),
            ("", "no FASTA records"),
            (">A\nAC\n>A\nGT\n", "duplicate"),
        ):
            with pytest.raises(SequenceError, match=msg):
                read_fasta_records(self._write(tmp_path, text, "bad.fa"))
            (tmp_path / "bad.fa").unlink()

    def test_lookup_verified_exact_length(self, tmp_path):
        p = self._write(tmp_path, f">TST1#LINE/TST\n{SEQ40}\n")
        res = library_lookup(p, "TST1", 40, library_name="my lib")
        assert res.status == STATUS_VERIFIED and res.verified
        assert res.sequence.sequence == SEQ40
        pv = res.provenance
        assert pv.provider == "user-provided library"
        assert pv.source_version == "my lib"
        assert pv.source_path == str(p)
        assert pv.checksum == sequence_checksum(SEQ40)
        assert pv.expected_model_length == 40

    def test_lookup_length_mismatch_incompatible_no_sequence(self, tmp_path):
        p = self._write(tmp_path, f">TST1\n{SEQ40}\n")
        res = library_lookup(p, "TST1", 41)
        assert res.status == STATUS_INCOMPATIBLE
        assert res.sequence is None  # never handed to display code
        assert "mismatch" in res.detail

    def test_lookup_missing_family_and_missing_file(self, tmp_path):
        p = self._write(tmp_path, ">A\nAC\n>B\nGT\n")
        assert library_lookup(p, "C", 2).status == STATUS_UNAVAILABLE
        res = library_lookup(tmp_path / "nope.fa", "A", 2)
        assert res.status == STATUS_UNAVAILABLE and "not found" in res.detail

    def test_lookup_single_record_file_used_for_any_family(self, tmp_path):
        p = self._write(tmp_path, f">whatever\n{SEQ40}\n")
        res = library_lookup(p, "TST1", 40, provider="user-provided FASTA")
        assert res.verified and res.provenance.matched_by == "only record in file"

    def test_lookup_invalid_alphabet_is_unavailable(self, tmp_path):
        p = self._write(tmp_path, ">TST1\nACGTX999\n")
        res = library_lookup(p, "TST1", 8)
        assert res.status == STATUS_UNAVAILABLE
        assert "invalid character" in res.detail


# --------------------------------------------------------------------- Dfam
class _FakeDfam:
    """Offline stand-in for DfamClient (same duck type)."""

    base_url = "https://dfam.test/api"

    def __init__(self, families, sequences, version="9.9"):
        self.families = families  # list of search results
        self.sequences = sequences  # accession -> fasta text
        self._version = version

    def version(self):
        return {"dfam_version": self._version, "release_date": "2099-01-01"}

    def search_families(self, name, limit=20):
        return [dict(f) for f in self.families]

    def family(self, accession):
        for f in self.families:
            if f["accession"] == accession:
                return dict(f)
        raise KeyError(accession)

    def family_sequence_fasta(self, accession):
        return self.sequences[accession]


class TestDfamSequence:
    def test_verified_fetch_with_provenance(self):
        client = _FakeDfam(
            [{"accession": "DF1", "name": "TST1", "length": 40}],
            {"DF1": f">DF1.1 TST1\n{SEQ40[:20]}\n{SEQ40[20:]}\n"},
        )
        res = sequence_lookup("TST1", 40, client=client)
        assert res.verified and res.sequence.sequence == SEQ40
        pv = res.provenance
        assert pv.accession == "DF1" and pv.matched_by == "name"
        assert pv.source_version == "Dfam 9.9"
        assert pv.checksum == sequence_checksum(SEQ40)

    def test_alias_match(self):
        client = _FakeDfam(
            [{"accession": "DF1", "name": "Other", "length": 40,
              "aliases": [{"alias": "TST1", "database": "Repbase"}]}],
            {"DF1": f">DF1\n{SEQ40}\n"},
        )
        res = sequence_lookup("TST1", 40, client=client)
        assert res.verified and res.provenance.matched_by == "alias:Repbase"

    def test_incompatible_records_every_rejection(self):
        # the L1HS shape: name-matching split fragments, neither compatible
        client = _FakeDfam(
            [{"accession": "DF1", "name": "TST1", "length": 9,
              "aliases": [{"alias": "TST1", "database": "Repbase"}]},
             {"accession": "DF2", "name": "TST1", "length": 21}],
            {},
        )
        res = sequence_lookup("TST1", 40, client=client)
        assert res.status == STATUS_INCOMPATIBLE and res.sequence is None
        assert "DF1" in res.detail and "DF2" in res.detail

    def test_no_candidates_unavailable(self):
        res = sequence_lookup("TST1", 40, client=_FakeDfam([], {}))
        assert res.status == STATUS_UNAVAILABLE
        assert "no Dfam family" in res.detail

    def test_record_length_lie_caught(self):
        # record metadata claims 40 bp but the fetched sequence is 39 bp:
        # the coordinate promise would be false -> incompatible
        client = _FakeDfam(
            [{"accession": "DF1", "name": "TST1", "length": 40}],
            {"DF1": f">DF1\n{SEQ40[:-1]}\n"},
        )
        res = sequence_lookup("TST1", 40, client=client)
        assert res.status == STATUS_INCOMPATIBLE and res.sequence is None
        assert "fetched sequence is 39 bp" in res.detail


# -------------------------------------------------------------------- cache
class TestSequenceCache:
    def _res(self):
        return SequenceResolution(
            status=STATUS_VERIFIED, detail="d",
            sequence=ConsensusSequence(SEQ40, SequenceProvenance(checksum="x")),
            provenance=SequenceProvenance(checksum="x"),
        )

    def test_roundtrip_and_drop(self, tmp_path):
        store_cached("TST1", 40, self._res(), cache_dir=tmp_path)
        back = load_cached("TST1", 40, cache_dir=tmp_path)
        assert back.verified and back.sequence.sequence == SEQ40
        assert drop_cached("TST1", 40, cache_dir=tmp_path)
        assert load_cached("TST1", 40, cache_dir=tmp_path) is None

    def test_identity_includes_model_length(self, tmp_path):
        # same family name, different coordinate model: distinct entries
        store_cached("TST1", 40, self._res(), cache_dir=tmp_path)
        assert load_cached("TST1", 41, cache_dir=tmp_path) is None
        assert cache_path("TST1", 40, cache_dir=tmp_path) != cache_path(
            "TST1", 41, cache_dir=tmp_path)

    def test_corrupt_cache_treated_as_missing(self, tmp_path):
        p = cache_path("TST1", 40, cache_dir=tmp_path)
        p.parent.mkdir(parents=True)
        p.write_text("{not json")
        assert load_cached("TST1", 40, cache_dir=tmp_path) is None


# ------------------------------------------------------------------ resolve
class TestResolve:
    def test_dfam_cache_first_offline(self, tmp_path):
        store_cached("TST1", 40, SequenceResolution(
            status=STATUS_VERIFIED, detail="cached snapshot",
            sequence=ConsensusSequence(SEQ40),
            provenance=SequenceProvenance()), cache_dir=tmp_path)
        res = resolve_sequence(SequenceConfig(enabled=True, source="dfam"),
                               "TST1", 40, cache_dir=tmp_path,
                               allow_network=False)
        assert res.verified and "(cached)" in res.detail

    def test_dfam_offline_uncached_unavailable(self, tmp_path):
        res = resolve_sequence(SequenceConfig(enabled=True, source="dfam"),
                               "TST1", 40, cache_dir=tmp_path,
                               allow_network=False)
        assert res.status == STATUS_UNAVAILABLE
        assert "not cached" in res.detail

    def test_dfam_network_error_degrades(self, tmp_path, monkeypatch):
        from tescope.core.features.dfam import DfamError
        from tescope.core.sequence import resolve as mod

        def boom(*a, **kw):
            raise DfamError("connection refused")

        monkeypatch.setattr(mod, "sequence_lookup", boom)
        res = resolve_sequence(SequenceConfig(enabled=True, source="dfam"),
                               "TST1", 40, cache_dir=tmp_path)
        assert res.status == STATUS_UNAVAILABLE
        assert "unreachable" in res.detail

    def test_dfam_lookup_result_is_cached(self, tmp_path, monkeypatch):
        from tescope.core.sequence import resolve as mod

        calls = []

        def fake(name, length, client=None):
            calls.append(name)
            return SequenceResolution(
                status=STATUS_VERIFIED, detail="live",
                sequence=ConsensusSequence(SEQ40),
                provenance=SequenceProvenance())

        monkeypatch.setattr(mod, "sequence_lookup", fake)
        cfg = SequenceConfig(enabled=True, source="dfam")
        r1 = resolve_sequence(cfg, "TST1", 40, cache_dir=tmp_path)
        r2 = resolve_sequence(cfg, "TST1", 40, cache_dir=tmp_path,
                              allow_network=False)
        assert r1.verified and r2.verified and calls == ["TST1"]

    def test_none_and_missing_path(self, tmp_path):
        assert resolve_sequence(SequenceConfig(), "T", 4).status == \
            STATUS_UNAVAILABLE
        res = resolve_sequence(
            SequenceConfig(enabled=True, source="library"), "T", 4)
        assert res.status == STATUS_UNAVAILABLE and "path" in res.detail

    def test_fasta_source(self, tmp_path):
        p = tmp_path / "one.fa"
        p.write_text(f">rec\n{SEQ40}\n")
        res = resolve_sequence(
            SequenceConfig(enabled=True, source="fasta", path=str(p)),
            "TST1", 40)
        assert res.verified
        assert res.provenance.provider == "user-provided FASTA"


# ---------------------------------------------------------------------- ops
class TestOps:
    SEQ = ConsensusSequence(
        SEQ40,
        SequenceProvenance(provider="Dfam", family="TST1", accession="DF1",
                           source_version="Dfam 9.9", checksum="abcd" * 16),
    )

    def test_subsequence_boundaries(self):
        assert subsequence(self.SEQ, 1, 40) == SEQ40
        assert subsequence(self.SEQ, 1, 1) == "A"
        assert subsequence(self.SEQ, 40, 40) == "A"
        for b1, b2 in ((0, 4), (1, 41), (5, 4)):
            with pytest.raises(RangeError):
                subsequence(self.SEQ, b1, b2)

    def test_fasta_partial_export_exact_bases(self):
        txt = fasta_text(self.SEQ, 9, 12)
        header, body = txt.strip().split("\n")
        assert body == SEQ40[8:12] == "GGGG"
        assert "positions=9-12(1-based)" in header
        assert "length=4" in header
        assert "Dfam:DF1" in header

    def test_fasta_full_export(self):
        txt = fasta_text(self.SEQ, tescope_version="9.9.9")
        assert "positions=1-40" in txt and "tescope=9.9.9" in txt
        assert "".join(txt.strip().split("\n")[1:]) == SEQ40

    def test_search_hits_1based(self):
        hits = search_literal(self.SEQ, "gggg")
        assert [(h.base_start, h.base_end) for h in hits] == [(9, 12)]
        assert SEQ40[hits[0].base_start - 1 : hits[0].base_end] == "GGGG"

    def test_search_overlapping_and_wildcard(self):
        hits = search_literal(ConsensusSequence("AAAA"), "AA")
        assert [(h.base_start) for h in hits] == [1, 2, 3]
        hits = search_literal(self.SEQ, "GNNC")
        assert (12, 15) in [(h.base_start, h.base_end) for h in hits]

    def test_search_revcomp_flagged(self):
        seq = ConsensusSequence("AACCGGTT")
        hits = search_literal(seq, "AACC", include_revcomp=True)
        # forward AACC at 1-4; reverse complement (GGTT) at 5-8, flagged
        assert [(h.base_start, h.base_end, h.strand) for h in hits] == \
            [(1, 4, "+"), (5, 8, "-")]
        assert "revcomp" in hits[1].label and "revcomp" not in hits[0].label
        # without the flag, only the forward hit is reported
        assert [(h.strand) for h in search_literal(seq, "AACC")] == ["+"]

    def test_search_errors(self):
        with pytest.raises(SequenceError):
            search_literal(self.SEQ, "")
        with pytest.raises(SequenceError):
            search_literal(self.SEQ, "AC-GT")
        assert search_literal(self.SEQ, "A" * 100) == []


# ------------------------------------------------------------- capabilities
class TestSequenceCapabilities:
    def test_none_resolution_all_disabled_with_reason(self):
        caps = sequence_capabilities(None)
        assert not caps.available and not caps.search
        assert "source" in caps.why_not("export_fasta")

    def test_verified_all_enabled(self):
        res = SequenceResolution(status=STATUS_VERIFIED,
                                 sequence=ConsensusSequence(SEQ40))
        caps = sequence_capabilities(res)
        assert caps.available and caps.base_display and caps.copy_sequence
        assert caps.export_fasta and caps.search and caps.verified

    def test_incompatible_disabled_with_model_reason(self):
        res = SequenceResolution(status=STATUS_INCOMPATIBLE,
                                 detail="length mismatch 40 vs 39")
        caps = sequence_capabilities(res)
        assert not caps.base_display
        why = caps.why_not("base_display")
        assert "does not match" in why and "length mismatch" in why


# ------------------------------------------------------------------- config
class TestSequenceConfig:
    def test_validation(self):
        SequenceConfig().validate()
        SequenceConfig(enabled=True, source="dfam").validate()
        with pytest.raises(ConfigError, match="requires sequence.path"):
            SequenceConfig(enabled=True, source="library").validate()
        with pytest.raises(ConfigError, match="sequence.source"):
            SequenceConfig(source="repbase-auto").validate()

    def test_analysis_config_roundtrip_and_backcompat(self):
        from tescope.config.model import RepeatSelection, TrackConfig

        cfg = AnalysisConfig(
            assembly="a", repeat=RepeatSelection(name="r"),
            tracks=[TrackConfig(file="x.bw")],
            sequence=SequenceConfig(enabled=True, source="fasta",
                                    path="lib.fa", library_name="mine"),
        )
        d = cfg.to_dict()
        back = AnalysisConfig.from_dict(d)
        assert back.sequence.enabled and back.sequence.path == "lib.fa"
        assert back.sequence.library_name == "mine"
        # configs written before the sequence section keep working
        d.pop("sequence")
        old = AnalysisConfig.from_dict(d)
        assert old.sequence.enabled is False
