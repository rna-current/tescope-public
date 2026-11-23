"""Locus table / subset / selection / estimate tests.

Expected values are HAND-CALCULATED from the golden dataset geometry
(tests/golden_dataset.py): the signal value of every base equals its
consensus position, so with bin_size=100 a full-length copy's bin j
mean is 100*j + 49.5 and its per-copy mean over observed bins is
499.5 — computed here on paper, not by re-running the implementation.
"""

from __future__ import annotations

import golden_dataset as gd
import numpy as np
import pytest

from tescope.config.model import AnalysisConfig, RepeatSelection, TrackConfig
from tescope.core.analysis import AnalysisError, run_analysis
from tescope.core.loci import (
    build_locus_table,
    export_selection,
    import_selection,
    load_estimates,
    locus_ids,
    subset_result,
)
from tescope.core.loci.estimates import EstimateError, match_estimates
from tescope.core.loci.subset import SubsetError


def make_config(golden_dir, **kw) -> AnalysisConfig:
    base = dict(
        assembly="testasm",
        repeat=RepeatSelection(name="TST1"),
        tracks=[TrackConfig(file=str(golden_dir["bigwig"]), name="sig"),
                TrackConfig(file=str(golden_dir["narrowpeak"]), name="peaks")],
    )
    base.update(kw)
    cfg = AnalysisConfig(**base)
    cfg.aggregation.bin_size = 100
    return cfg


@pytest.fixture(scope="module")
def result(golden_dir, golden_store):
    return run_analysis(make_config(golden_dir), store=golden_store)


# -- locus table -------------------------------------------------------------


def test_locus_ids_are_stable_coordinate_strings(result):
    ids = locus_ids(result.copies)
    gs, ge, strand = gd.COPIES["F"][0], gd.COPIES["F"][1], gd.COPIES["F"][2]
    assert f"{gd.CHROM}:{gs}-{ge}({strand})" in set(ids)
    assert len(ids) == len(set(ids))  # unique for this dataset


def test_locus_table_hand_calculated_means(result):
    table = build_locus_table(result)
    by_id = table.set_index("locus_id")
    # copy F: full length, 10 observed 100-bp bins, mean_j = 100j+49.5
    f_id = f"{gd.CHROM}:10000-11000(+)"
    assert by_id.loc[f_id, "obs_bins_sig"] == 10
    assert by_id.loc[f_id, "mean_sig"] == pytest.approx(499.5)
    # copy R: reverse strand maps to the identical consensus profile
    assert by_id.loc[f"{gd.CHROM}:20000-21000(-)", "mean_sig"] == pytest.approx(499.5)
    # copy SHORT: cons 100-150 -> one observed bin, values 100..149
    assert by_id.loc[f"{gd.CHROM}:60000-60050(+)", "obs_bins_sig"] == 1
    assert by_id.loc[f"{gd.CHROM}:60000-60050(+)", "mean_sig"] == pytest.approx(124.5)
    # copy MID: bins 4 and 5 -> (449.5 + 549.5)/2
    assert by_id.loc[f"{gd.CHROM}:50000-50200(+)", "mean_sig"] == pytest.approx(499.5)


def test_locus_table_carries_identity_and_classes(result):
    table = build_locus_table(result)
    assert set(
        ["locus_id", "chrom", "genome_start", "genome_end", "strand",
         "completeness", "consensus_coverage", "rm_id", "rm_id_fragments"]
    ) <= set(table.columns)
    by_id = table.set_index("locus_id")
    for name, (gs, ge, strand, *_rest) in gd.COPIES.items():
        lid = f"{gd.CHROM}:{gs}-{ge}({strand})"
        if lid in by_id.index:  # STAR may be filtered by include_secondary
            assert by_id.loc[lid, "completeness"] == gd.EXPECTED_COMPLETENESS[name]


def test_locus_table_peak_counts_are_observed_zeros_not_missing(result):
    table = build_locus_table(result)
    col = table["peaks_peaks"]
    # the track covers chrT, so every value is a real number (0 included)
    assert col.notna().all()
    assert (col >= 0).all()


# -- subset analysis ---------------------------------------------------------


def test_subset_hand_calculated_and_original_untouched(result):
    before = result.continuous_tracks()[0].profile["mean"].copy()
    ids = [f"{gd.CHROM}:10000-11000(+)", f"{gd.CHROM}:20000-21000(-)"]
    sub = subset_result(result, ids)
    prof = sub.profiles["sig"]
    # F and R have identical consensus rows: mean_j = 100j+49.5, n=2
    assert prof["mean"].to_numpy() == pytest.approx(
        np.arange(10) * 100 + 49.5
    )
    assert (prof["n_observations"] == 2).all()
    assert (prof["n_elements"] == 2).all()
    assert sub.n_subset == 2 and sub.n_eligible == len(result.copies)
    # the eligible population's aggregate is untouched
    after = result.continuous_tracks()[0].profile["mean"]
    assert after.to_numpy() == pytest.approx(before.to_numpy(), nan_ok=True)


def test_subset_reports_invalid_ids_and_requires_two(result):
    good = f"{gd.CHROM}:10000-11000(+)"
    sub = subset_result(result, [good, f"{gd.CHROM}:20000-21000(-)",
                                 "chrZ:1-2(+)"])
    assert sub.invalid_ids == ["chrZ:1-2(+)"]
    with pytest.raises(SubsetError):
        subset_result(result, [good, "chrZ:1-2(+)"])


def test_subset_occurrence_uses_subset_denominator(result):
    ids = [f"{gd.CHROM}:10000-11000(+)", f"{gd.CHROM}:20000-21000(-)"]
    sub = subset_result(result, ids)
    occ = sub.occurrences["peaks"]
    assert (occ["n_elements"] <= 2).all()


# -- selection files ---------------------------------------------------------


def test_selection_round_trip(result, tmp_path):
    ids = [f"{gd.CHROM}:10000-11000(+)", f"{gd.CHROM}:30000-30400(+)"]
    p = tmp_path / "sel.yaml"
    export_selection(result, ids, str(p))
    valid, invalid, warn = import_selection(str(p), result)
    assert valid == ids and invalid == [] and warn == ""


def test_selection_import_reports_unknown_ids(result, tmp_path):
    p = tmp_path / "sel.yaml"
    export_selection(result, [f"{gd.CHROM}:10000-11000(+)"], str(p))
    text = p.read_text().replace("10000-11000", "10001-11001")
    p.write_text(text)
    valid, invalid, _ = import_selection(str(p), result)
    assert valid == [] and invalid == [f"{gd.CHROM}:10001-11001(+)"]


def test_selection_assembly_mismatch_fails(result, tmp_path):
    p = tmp_path / "sel.yaml"
    export_selection(result, [f"{gd.CHROM}:10000-11000(+)"], str(p))
    p.write_text(p.read_text().replace("assembly: testasm", "assembly: hg38"))
    from tescope.core.loci import SelectionError

    with pytest.raises(SelectionError):
        import_selection(str(p), result)


# -- upstream estimates ------------------------------------------------------


def _write_estimates(tmp_path, body, header_extra=""):
    p = tmp_path / "est.tsv"
    p.write_text(
        "# tescope-locus-estimates v1\n"
        "# tool: telescope 1.0.3\n"
        "# assay: RNA-seq\n"
        "# units: fragments\n"
        "# total_includes_unique: true\n"
        f"{header_extra}"
        "locus_id\ttotal\tunique\tscore\n" + body
    )
    return str(p)


def test_estimates_load_and_match(result, tmp_path):
    f_id = f"{gd.CHROM}:10000-11000(+)"
    path = _write_estimates(
        tmp_path,
        f"{f_id}\t120\t45\t0.97\nchrZ:5-9(+)\t3\t1\t0.5\n",
    )
    est = load_estimates(path)
    assert est.total_includes_unique is True and est.assay == "RNA-seq"
    table = build_locus_table(result)
    merged, report = match_estimates(est, table, "testasm")
    assert report.n_matched == 1 and report.n_unmatched == 1
    assert report.unmatched_ids == ["chrZ:5-9(+)"]
    row = merged.set_index("locus_id").loc[f_id]
    assert row["est_total"] == 120 and row["est_unique"] == 45
    # unmatched analyzed loci: unavailable (NaN), never zero
    other = merged.set_index("locus_id").loc[f"{gd.CHROM}:20000-21000(-)"]
    assert np.isnan(other["est_total"])
    # the NEVER-ADD invariant: no derived total+unique column exists
    assert not any("sum" in c or "combined" in c for c in merged.columns)


def test_estimates_reject_inconsistent_decomposition(result, tmp_path):
    path = _write_estimates(tmp_path, f"{gd.CHROM}:10000-11000(+)\t10\t45\t0.9\n")
    with pytest.raises(EstimateError, match="unique > total"):
        load_estimates(path)


def test_estimates_require_semantics_directives(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text(
        "# tescope-locus-estimates v1\n# assay: RNA-seq\n"
        "locus_id\ttotal\nchr1:1-2(+)\t5\n"
    )
    with pytest.raises(EstimateError, match="units"):
        load_estimates(str(p))


def test_estimates_assembly_mismatch_fails(result, tmp_path):
    path = _write_estimates(
        tmp_path, f"{gd.CHROM}:10000-11000(+)\t120\t45\t0.9\n",
        header_extra="# assembly: hg38\n",
    )
    est = load_estimates(path)
    with pytest.raises(EstimateError, match="assembly"):
        match_estimates(est, build_locus_table(result), "testasm")


def test_estimates_through_run_analysis(golden_dir, golden_store, tmp_path):
    f_id = f"{gd.CHROM}:10000-11000(+)"
    path = _write_estimates(tmp_path, f"{f_id}\t120\t45\t0.97\n")
    cfg = make_config(golden_dir)
    cfg.loci.estimates_file = path
    res = run_analysis(cfg, store=golden_store)
    table = res.locus_table()
    assert "est_total" in table.columns
    assert res.estimate_report.n_matched == 1
    # positional outputs are untouched by estimates
    assert "est_total" not in res.continuous_tracks()[0].profile.columns


def test_run_analysis_fails_loudly_on_bad_estimate_file(
    golden_dir, golden_store, tmp_path
):
    bad = tmp_path / "bad.tsv"
    bad.write_text("not an estimates file\n")
    cfg = make_config(golden_dir)
    cfg.loci.estimates_file = str(bad)
    with pytest.raises(AnalysisError, match="estimate"):
        run_analysis(cfg, store=golden_store)
