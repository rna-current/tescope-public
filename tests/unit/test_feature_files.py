"""Feature-file import (TSV template / BED / GFF3) + locus projection."""

from __future__ import annotations

import golden_dataset as gd
import pytest

from tescope.config.model import AnalysisConfig, RepeatSelection, TrackConfig
from tescope.core.analysis import AnalysisError, run_analysis
from tescope.core.coordinates.mapping import (
    CopyAlignment,
    consensus_interval_to_genomic,
    genomic_interval_to_consensus,
)
from tescope.core.features.files import (
    FeatureFileError,
    attach_genome_features,
    parse_feature_file,
)

TEMPLATE = """# tescope-features v1
# target: consensus
# reference: TST1
# coordinates: 0-based-half-open
feature_id\tstart\tend\tstrand\tlabel\ttype\tparent_id\tsource
ORF1\t100\t400\t+\tORF1\tCDS\t\tliterature
UTR5\t0\t100\t+\t5'UTR\tUTR\t\tliterature
SUB\t150\t200\t+\tsubdomain\tdomain\tORF1\tliterature
"""


def _cfg(golden_dir, **kw):
    cfg = AnalysisConfig(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file=str(golden_dir["bigwig"]), name="sig")],
        **kw,
    )
    cfg.aggregation.bin_size = 100
    return cfg


# -- template TSV ------------------------------------------------------------


def test_template_consensus_parse(tmp_path):
    p = tmp_path / "f.tsv"
    p.write_text(TEMPLATE)
    r = parse_feature_file(str(p))
    assert r.target == "consensus" and r.reference == "TST1"
    assert [f.feature_id for f in r.features] == ["ORF1", "UTR5", "SUB"]
    assert r.features[2].parent_id == "ORF1"


def test_template_rejects_amino_acid_coordinates(tmp_path):
    p = tmp_path / "f.tsv"
    p.write_text(TEMPLATE.replace("# target:", "# units: aa\n# target:"))
    with pytest.raises(FeatureFileError, match="amino-acid"):
        parse_feature_file(str(p))


def test_template_validation_collects_errors(tmp_path):
    p = tmp_path / "f.tsv"
    p.write_text(
        "# tescope-features v1\n# target: consensus\n"
        "feature_id\tstart\tend\tstrand\tparent_id\n"
        "A\t100\t50\t+\t\n"  # end <= start
        "A\t10\t20\t?\tGHOST\n"  # dup id, bad strand, dangling parent
    )
    with pytest.raises(FeatureFileError) as exc:
        parse_feature_file(str(p))
    msg = str(exc.value)
    for frag in ("end (50) must be > start", "duplicate feature_id",
                 "strand", "parent_id 'GHOST'"):
        assert frag in msg


def test_template_genome_target_requires_contig(tmp_path):
    p = tmp_path / "f.tsv"
    p.write_text(
        "# tescope-features v1\n# target: genome\n# reference: testasm\n"
        "feature_id\tstart\tend\n"
        f"{gd.CHROM}:geneA\t9000\t9500\n"
    )
    r = parse_feature_file(str(p))
    assert r.features[0].chrom == gd.CHROM
    assert r.features[0].feature_id == "geneA"
    p.write_text(
        "# tescope-features v1\n# target: genome\n# reference: testasm\n"
        "feature_id\tstart\tend\nnocontig\t1\t2\n"
    )
    with pytest.raises(FeatureFileError, match="contig"):
        parse_feature_file(str(p))


# -- BED / GFF3 conventions --------------------------------------------------


def test_bed_is_zero_based_half_open_by_spec(tmp_path):
    p = tmp_path / "f.bed"
    p.write_text(f"{gd.CHROM}\t100\t200\tgeneA\t0\t+\n")
    r = parse_feature_file(str(p))
    f = r.features[0]
    assert (f.start, f.end, f.strand, f.chrom) == (100, 200, "+", gd.CHROM)
    assert "BED specification" in r.source_convention


def test_gff3_one_based_inclusive_converted(tmp_path):
    p = tmp_path / "f.gff3"
    p.write_text(
        "##gff-version 3\n"
        f"{gd.CHROM}\tsrc\tgene\t101\t200\t.\t+\t.\tID=geneA;Name=GeneA\n"
        f"{gd.CHROM}\tsrc\texon\t101\t150\t.\t+\t.\tID=ex1;Parent=geneA\n"
    )
    r = parse_feature_file(str(p))
    f = r.features[0]
    # GFF3 101..200 (1-based inclusive) == internal [100, 200)
    assert (f.start, f.end) == (100, 200)
    assert r.features[1].parent_id == "geneA"
    assert "GFF3 specification" in r.source_convention


def test_bed_duplicate_names_suffixed_not_merged(tmp_path):
    p = tmp_path / "f.bed"
    p.write_text(
        f"{gd.CHROM}\t10\t20\tX\t0\t+\n{gd.CHROM}\t30\t40\tX\t0\t+\n"
    )
    r = parse_feature_file(str(p))
    assert [f.feature_id for f in r.features] == ["X", "X#2"]
    assert any("NOT merged" in w for w in r.warnings)


# -- genome attachment gates -------------------------------------------------


def test_genome_attach_assembly_and_contig_gates(tmp_path):
    p = tmp_path / "f.bed"
    p.write_text(f"{gd.CHROM}\t100\t200\tgeneA\n")
    parsed = parse_feature_file(str(p))
    with pytest.raises(FeatureFileError, match="assembly"):
        attach_genome_features(parsed, "testasm", "hg38", [gd.CHROM])
    with pytest.raises(FeatureFileError, match="no assembly declared"):
        attach_genome_features(parsed, "testasm", "", [gd.CHROM])
    p2 = tmp_path / "g.bed"
    p2.write_text("chrUnknown\t1\t2\tz\n")
    parsed2 = parse_feature_file(str(p2))
    with pytest.raises(FeatureFileError, match="contig"):
        attach_genome_features(parsed2, "testasm", "testasm", [gd.CHROM])
    gset = attach_genome_features(parsed, "testasm", "testasm", [gd.CHROM])
    assert gset.overlapping(gd.CHROM, 150, 160)
    assert not gset.overlapping(gd.CHROM, 500, 600)
    assert gset.overlapping(gd.CHROM, 250, 300, flank=100)


# -- run_analysis integration ------------------------------------------------


def test_feature_files_through_run_analysis(golden_dir, golden_store, tmp_path):
    cons = tmp_path / "cons.tsv"
    cons.write_text(TEMPLATE)
    gen = tmp_path / "genes.bed"
    gen.write_text(f"{gd.CHROM}\t9000\t9500\tneighborGene\t0\t+\n")
    cfg = _cfg(golden_dir)
    cfg.annotations.feature_files = [
        {"path": str(cons)},
        {"path": str(gen), "reference": "testasm"},
    ]
    res = run_analysis(cfg, store=golden_store)
    assert res.features is not None
    assert {f.name for f in res.features.user_features} == {
        "ORF1", "5'UTR", "subdomain"
    }
    assert len(res.genome_features) == 1
    assert res.genome_features[0].features[0].feature_id == "neighborGene"


def test_wrong_consensus_model_withheld(golden_dir, golden_store, tmp_path):
    cons = tmp_path / "cons.tsv"
    cons.write_text(TEMPLATE.replace("# reference: TST1", "# reference: L1HS"))
    cfg = _cfg(golden_dir)
    cfg.annotations.feature_files = [{"path": str(cons)}]
    with pytest.raises(AnalysisError, match="L1HS"):
        run_analysis(cfg, store=golden_store)


def test_out_of_bounds_consensus_feature_fails(golden_dir, golden_store, tmp_path):
    cons = tmp_path / "cons.tsv"
    cons.write_text(
        "# tescope-features v1\n# target: consensus\n# reference: TST1\n"
        "feature_id\tstart\tend\nTOObig\t0\t5000\n"
    )
    cfg = _cfg(golden_dir)
    cfg.annotations.feature_files = [{"path": str(cons)}]
    with pytest.raises(AnalysisError):
        run_analysis(cfg, store=golden_store)


# -- consensus -> locus projection (inverse mapping) -------------------------


def test_inverse_mapping_round_trip_and_strand():
    # indel-free copy, + strand: exact inverse of the forward mapping
    c = CopyAlignment(10_000, 11_000, 0, 1000, "+")
    fwd = genomic_interval_to_consensus(c, 10_100, 10_400)
    assert fwd == (100, 400)
    inv = consensus_interval_to_genomic(c, 100, 400)
    assert inv == (10_100, 10_400, False)
    # reverse strand: consensus start maps to the genomic END
    r = CopyAlignment(20_000, 21_000, 0, 1000, "-")
    assert consensus_interval_to_genomic(r, 0, 100) == (20_900, 21_000, False)
    assert consensus_interval_to_genomic(r, 900, 1000) == (20_000, 20_100, False)


def test_inverse_mapping_truncation_and_outside():
    # P3R-like: 5'-truncated reverse copy covering consensus 700-1000
    c = CopyAlignment(40_000, 40_300, 700, 1000, "-")
    # a feature fully 5' of the footprint has NO verified correspondence
    assert consensus_interval_to_genomic(c, 100, 400) is None
    # a feature straddling the edge projects only its overlap, flagged
    g0, g1, partial = consensus_interval_to_genomic(c, 650, 750)
    assert partial is True
    assert (g0, g1) == (40_250, 40_300)  # cons 700-750 -> genomic tail


def test_inverse_mapping_one_bp_feature_with_indels():
    # 2:1 scale (genomic length 200, consensus span 100)
    c = CopyAlignment(1_000, 1_200, 0, 100, "+")
    g0, g1, partial = consensus_interval_to_genomic(c, 50, 51)
    assert not partial
    assert g0 == 1_100 and g1 == 1_102  # one consensus base ~ two genomic
