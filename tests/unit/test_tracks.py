"""Track readers: bigWig (all backends), bedGraph, BED/narrowPeak, points, loader."""

from __future__ import annotations

import golden_dataset as gd
import numpy as np
import pytest

from tescope.core.chromnames import ChromMismatchError, reconcile, require_overlap
from tescope.core.tracks.base import TrackKind, TrackMeta
from tescope.core.tracks.bed import open_intervals, open_points
from tescope.core.tracks.bedgraph import BedGraphParseError, open_bedgraph
from tescope.core.tracks.bigwig import BigWigTrack, available_backends
from tescope.core.tracks.loader import UnknownTrackFormatError, detect_format, open_track


class TestBedGraph:
    def test_values(self, golden_dir):
        t = open_bedgraph(golden_dir["bedgraph"])
        v = t.values(gd.CHROM, 10_000, 10_010)
        assert np.allclose(v, np.arange(10))

    def test_nan_outside_data(self, golden_dir):
        t = open_bedgraph(golden_dir["bedgraph"])
        v = t.values(gd.CHROM, 5_000, 5_010)
        assert np.isnan(v).all()

    def test_missing_chrom_is_nan(self, golden_dir):
        t = open_bedgraph(golden_dir["bedgraph"])
        assert np.isnan(t.values("chrOther", 0, 5)).all()

    def test_rejects_overlaps(self, tmp_path):
        p = tmp_path / "bad.bedGraph"
        p.write_text("chr1\t0\t10\t1.0\nchr1\t5\t15\t2.0\n")
        with pytest.raises(BedGraphParseError, match="overlap"):
            open_bedgraph(str(p))

    def test_rejects_malformed(self, tmp_path):
        p = tmp_path / "bad2.bedGraph"
        p.write_text("chr1\t0\tx\t1.0\n")
        with pytest.raises(BedGraphParseError, match="bad2"):
            open_bedgraph(str(p))


class TestBigWig:
    def test_backends_available(self):
        assert available_backends(), "no bigWig backend importable"

    @pytest.mark.parametrize("backend", available_backends())
    def test_values_per_backend(self, golden_dir, backend):
        meta = TrackMeta(name="s", path=golden_dir["bigwig"], kind=TrackKind.CONTINUOUS,
                         format="bigwig")
        t = BigWigTrack(meta, backend=backend)
        v = t.values(gd.CHROM, 10_000, 10_010)
        assert np.allclose(v, np.arange(10)), backend
        assert np.isnan(t.values(gd.CHROM, 5_000, 5_005)).all()
        t.close()

    def test_backend_parity_with_bedgraph(self, golden_dir):
        bg = open_bedgraph(golden_dir["bedgraph"])
        expected = bg.values(gd.CHROM, 9_990, 11_010)
        for backend in available_backends():
            meta = TrackMeta(name="s", path=golden_dir["bigwig"], kind=TrackKind.CONTINUOUS,
                             format="bigwig")
            t = BigWigTrack(meta, backend=backend)
            got = t.values(gd.CHROM, 9_990, 11_010)
            np.testing.assert_allclose(got, expected, equal_nan=True, err_msg=backend)
            t.close()


class TestIntervals:
    def test_fetch(self, golden_dir):
        t = open_intervals(golden_dir["peaks_bed"])
        ivs = t.fetch(gd.CHROM, 10_000, 11_000)
        assert [iv.name for iv in ivs] == ["peak_F"]
        assert (ivs[0].start, ivs[0].end) == gd.PEAKS["peak_F"]
        assert t.count() == len(gd.PEAKS)

    def test_narrowpeak_signalvalue_as_score(self, golden_dir):
        t = open_intervals(golden_dir["narrowpeak"])
        ivs = t.fetch(gd.CHROM, 10_000, 11_000)
        assert ivs[0].score == pytest.approx(7.5)

    def test_overlap_semantics(self, golden_dir):
        t = open_intervals(golden_dir["peaks_bed"])
        # query window strictly inside the peak still finds it
        ivs = t.fetch(gd.CHROM, 10_260, 10_270)
        assert len(ivs) == 1
        assert t.fetch(gd.CHROM, 10_300, 10_400) == []  # half-open: no touch


class TestPoints:
    def test_fetch(self, golden_dir):
        t = open_points(golden_dir["sites"])
        pos, vals = t.fetch(gd.CHROM, 0, gd.CHROM_LEN)
        assert pos.tolist() == sorted(gd.SITES)
        assert t.count() == len(gd.SITES)
        assert vals.tolist() == [1.0, 2.0]


class TestLoader:
    @pytest.mark.parametrize(
        "fname,kind,fmt",
        [
            ("x.bw", "continuous", "bigwig"),
            ("x.bigWig", "continuous", "bigwig"),
            ("x.bedGraph", "continuous", "bedgraph"),
            ("x.bdg.gz", "continuous", "bedgraph"),
            ("x.bed", "interval", "bed"),
            ("x.narrowPeak", "interval", "narrowpeak"),
            ("x.broadPeak.gz", "interval", "broadpeak"),
            ("x.bam", "continuous", "bam"),
        ],
    )
    def test_detect(self, fname, kind, fmt):
        k, f = detect_format(fname)
        assert (k.value, f) == (kind, fmt)

    def test_unknown(self):
        with pytest.raises(UnknownTrackFormatError):
            detect_format("mystery.xyz")

    def test_open_as_points(self, golden_dir):
        t = open_track(golden_dir["sites"], kind="point")
        assert t.meta.kind is TrackKind.POINT


class TestChromNames:
    def test_exact_match(self):
        r = reconcile(["chr1", "chr2"], ["chr1", "chr2", "chr3"])
        assert r.matched == {"chr1": "chr1", "chr2": "chr2"}
        assert r.translate("chr1") == "chr1"
        assert r.translate("chr3") is None

    def test_alias(self):
        r = reconcile(["1", "MT"], ["chr1", "chrM"])
        assert r.matched == {"1": "chr1", "MT": "chrM"}
        assert r.aliased == r.matched
        assert r.translate("chrM") == "MT"

    def test_alias_disabled(self):
        r = reconcile(["1"], ["chr1"], allow_alias=False)
        assert r.unmatched_track == ["1"]

    def test_no_silent_remap_of_unknown(self):
        r = reconcile(["scaffold_1"], ["chr1"])
        assert r.unmatched_track == ["scaffold_1"]

    def test_require_overlap_raises(self):
        r = reconcile(["weird"], ["chr1"])
        with pytest.raises(ChromMismatchError, match="weird"):
            require_overlap(r, "track.bw")
