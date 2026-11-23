"""BAM coverage track (optional pysam path): filters, index handling."""

from __future__ import annotations

import pytest

pysam = pytest.importorskip("pysam")

import golden_dataset as gd  # noqa: E402

from tescope.core.tracks.bam import BamCoverageError, open_bam_coverage  # noqa: E402

CHROM = gd.CHROM


def _make_bam(path, reads):
    """reads: list of (start, length, mapq, is_reverse, is_duplicate)."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": gd.CHROM_LEN}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, (start, length, mapq, rev, dup) in enumerate(sorted(reads)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * length
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = mapq
            a.cigarstring = f"{length}M"
            a.flag = (0x10 if rev else 0) | (0x400 if dup else 0)
            a.query_qualities = pysam.qualitystring_to_array("I" * length)
            bam.write(a)
    pysam.index(str(path))


@pytest.fixture(scope="module")
def bam_path(tmp_path_factory):
    p = tmp_path_factory.mktemp("bam") / "reads.bam"
    _make_bam(
        p,
        [
            (100, 50, 60, False, False),   # forward, high mapq
            (120, 50, 60, True, False),    # reverse, high mapq
            (100, 50, 0, False, False),    # mapq 0 (multi-mapper-like)
            (100, 50, 60, False, True),    # duplicate
        ],
    )
    return str(p)


class TestBamCoverage:
    def test_basic_coverage(self, bam_path):
        t = open_bam_coverage(bam_path, min_mapq=1)
        v = t.values(CHROM, 100, 170)
        # dup excluded by default, mapq0 excluded: fwd covers 100-150, rev 120-170
        assert v[0] == 1.0            # only fwd
        assert v[25] == 2.0           # fwd + rev overlap 120-150
        assert v[55] == 1.0           # only rev
        t.close()

    def test_mapq_zero_included_when_allowed(self, bam_path):
        t = open_bam_coverage(bam_path, min_mapq=0)
        v = t.values(CHROM, 100, 120)
        assert v[0] == 2.0            # fwd + mapq0
        t.close()

    def test_duplicates_optional(self, bam_path):
        t = open_bam_coverage(bam_path, min_mapq=1, exclude_duplicates=False)
        assert t.values(CHROM, 100, 101)[0] == 2.0  # fwd + dup
        t.close()

    def test_strand_selection(self, bam_path):
        t = open_bam_coverage(bam_path, min_mapq=1, strand="-")
        v = t.values(CHROM, 100, 170)
        assert v[0] == 0.0 and v[25] == 1.0
        t.close()

    def test_missing_index_message(self, tmp_path, bam_path):
        import shutil

        p = tmp_path / "noindex.bam"
        shutil.copy(bam_path, p)
        with pytest.raises(BamCoverageError, match="index missing"):
            open_bam_coverage(str(p))

    def test_auto_index(self, tmp_path, bam_path):
        import shutil

        p = tmp_path / "autoindex.bam"
        shutil.copy(bam_path, p)
        t = open_bam_coverage(str(p), create_index=True, min_mapq=1)
        assert t.values(CHROM, 100, 101)[0] == 1.0
        t.close()

    def test_caveat_recorded(self, bam_path):
        t = open_bam_coverage(bam_path)
        assert "multi-mapping" in t.meta.metadata["caveat"]
        t.close()
