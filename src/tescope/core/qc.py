"""QC summary generation: one tidy table per analysis, GUI-showable and exportable."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .analysis import AnalysisResult
from .tracks.base import TrackKind


def qc_table(result: AnalysisResult) -> pd.DataFrame:
    """Tidy QC table with columns: section, metric, value."""
    rows: list[tuple[str, str, str]] = []

    def add(section: str, metric: str, value) -> None:
        rows.append((section, metric, str(value)))

    add("analysis", "assembly", result.config.assembly)
    add("analysis", "repeat_name", result.repeat_name)
    add("analysis", "consensus_length_bp", result.cons_length)
    add("analysis", "bin_size_bp", result.config.aggregation.bin_size)
    add("analysis", "n_bins", len(result.bin_edges) - 1)
    prov = result.provenance
    if prov:
        add("annotation", "source", prov.get("source", ""))
        add("annotation", "source_url_or_path", prov.get("source_url") or prov.get("source_path", ""))
        add("annotation", "retrieved_at", prov.get("retrieved_at", ""))
        add("annotation", "md5", prov.get("md5", ""))

    for label, n in result.filter_report.to_rows():
        add("filtering", label, n)
    add("filtering", "Filtered out", result.n_annotated - result.n_selected)

    if result.n_selected:
        comp = result.copies["completeness"].value_counts()
        for cls, n in comp.items():
            add("completeness", str(cls), int(n))
        strands = result.copies["strand"].value_counts()
        for s, n in strands.items():
            add("strand", f"copies on {s}", int(n))

    n_elements = result.element_covered.sum(axis=0)
    if len(n_elements):
        add("consensus coverage", "max elements per bin", int(n_elements.max()))
        add("consensus coverage", "min elements per bin", int(n_elements.min()))
        add("consensus coverage", "median elements per bin", float(np.median(n_elements)))
        frac_covered = float((n_elements > 0).mean())
        add("consensus coverage", "fraction of bins covered by >=1 element", f"{frac_covered:.3f}")

    add("tracks", "n_tracks", len(result.tracks))
    for t in result.tracks:
        sec = f"track: {t.name}"
        add(sec, "kind", t.kind.value)
        add(sec, "file", t.config.file)
        if t.config.group:
            add(sec, "group", t.config.group)
        if t.config.replicate is not None:
            add(sec, "replicate", t.config.replicate)
        if t.chrom_report:
            add(sec, "chromosomes", t.chrom_report.summary())
        if t.kind is TrackKind.CONTINUOUS:
            is_bam = (t.config.format == "bam"
                      or t.config.file.lower().endswith((".bam", ".cram")))
            if is_bam:
                mc = t.config.mapping_class or "all"
                pol = f"class={mc}, method={t.config.mapping_method}"
                if t.config.mapping_method == "mapq" and mc != "all":
                    pol += f" (unique_mapq={t.config.unique_mapq})"
                if t.config.min_mapq:
                    pol += f", min_mapq={t.config.min_mapq}"
                add(sec, "mapping policy", pol)
                if mc == "unique":
                    add(sec, "mapping note",
                        "unique-only coverage: multi-mapping reads excluded "
                        "by explicit policy (aligner-dependent classification)")
            else:
                add(sec, "mapping policy", "unavailable from track "
                    "(pre-aggregated signal; no read-level information)")
        if t.kind is TrackKind.CONTINUOUS and t.matrix is not None:
            obs = np.isfinite(t.matrix.values)
            covered = t.matrix.covered
            with np.errstate(invalid="ignore"):
                missing = 1.0 - (obs.sum() / covered.sum()) if covered.sum() else float("nan")
            add(sec, "copies skipped (chromosome absent)", t.matrix.n_skipped_chrom)
            add(sec, "signal missingness (covered bins without data)", f"{missing:.3f}")
            if t.config.format == "bam" or t.config.file.lower().endswith((".bam", ".cram")):
                add(sec, "caveat", "BAM coverage: multi-mapping reads in TE regions not rescued")
        if t.kind is TrackKind.INTERVAL:
            add(sec, "peaks in file", t.n_source_items)
            add(sec, "projected (copy, peak) pairs", 0 if t.projected is None else len(t.projected))
            if t.projected is not None and len(t.projected):
                add(sec, "TE copies with >=1 peak", t.projected["copy_idx"].nunique())
        if t.kind is TrackKind.POINT:
            add(sec, "points in file", t.n_source_items)
            add(sec, "points inside TE copies", 0 if t.points is None else len(t.points))

    for g, prof in result.group_profiles.items():
        add("groups", g, f"{int(prof['n_replicates'].max())} replicate track(s)")

    if result.features is not None:
        fs = result.features
        add("annotation", "status", fs.status)
        if fs.provenance is not None:
            pv = fs.provenance
            add("annotation", "source", f"{pv.source} {pv.database_version}".strip())
            if pv.accession:
                add("annotation", "accession", f"{pv.accession} ({pv.record_name})")
            add("annotation", "compatibility", pv.compatibility)
            if pv.compatibility_detail:
                add("annotation", "compatibility detail", pv.compatibility_detail)
            if pv.retrieved_at:
                add("annotation", "retrieved", pv.retrieved_at)
        add("annotation", "database features", len(fs.database_features))
        add("annotation", "user-defined features", len(fs.user_features))
        # Feature coverage QC — definitions (docs/annotations.md):
        #   overlapping    = copy cons interval intersects the feature interval
        #   fully covering = copy cons_start <= f.start and cons_end >= f.end
        # Denominator = analyzed (post-filter) copies. Purely descriptive;
        # does not alter any profile calculation.
        n = max(result.n_selected, 1)
        cs = result.copies["cons_start"].to_numpy()
        ce = result.copies["cons_end"].to_numpy()
        for f in fs.all_features():
            over = int(((cs < f.end) & (ce > f.start)).sum())
            full = int(((cs <= f.start) & (ce >= f.end)).sum())
            sec = f"feature: {f.name}"
            add(sec, "origin", f.origin)
            add(sec, "consensus interval", f"{f.start}-{f.end}")
            add(sec, "copies overlapping", f"{over} ({over / n:.1%} of analyzed)")
            add(sec, "copies fully covering", f"{full} ({full / n:.1%} of analyzed)")

    if getattr(result, "sequence", None) is not None:
        sq = result.sequence
        add("sequence", "status", sq.status)
        if sq.detail:
            add("sequence", "detail", sq.detail)
        if sq.provenance is not None:
            pv = sq.provenance
            add("sequence", "provider", pv.provider)
            if pv.accession:
                add("sequence", "accession", f"{pv.accession} ({pv.record_name})")
            elif pv.record_name:
                add("sequence", "record", pv.record_name)
            if pv.source_version:
                add("sequence", "source release", pv.source_version)
            if pv.source_path:
                add("sequence", "source path", pv.source_path)
            if pv.sequence_length is not None:
                add("sequence", "sequence length",
                    f"{pv.sequence_length:,} bp (model "
                    f"{pv.expected_model_length or '?'} bp)")
            if pv.checksum:
                add("sequence", "sha256", pv.checksum)
            if pv.retrieved_at:
                add("sequence", "retrieved", pv.retrieved_at)
            add("sequence", "coordinate convention",
                "displayed base positions are 1-based inclusive")

    for gset in result.genome_features:
        sec = "genome features"
        add(sec, "file", gset.path)
        add(sec, "assembly", gset.assembly)
        add(sec, "features", f"{len(gset.features):,} ({gset.format}, "
            "stored 0-based half-open)")

    if result.locus_estimates is not None:
        est = result.locus_estimates
        sec = "locus estimates"
        add(sec, "file", est.path)
        add(sec, "tool", est.tool or "unknown (not declared)")
        add(sec, "assay", est.assay)
        add(sec, "units", est.units)
        tiu = est.total_includes_unique
        add(sec, "total includes unique",
            "unknown (declared)" if tiu is None else str(tiu).lower())
        add(sec, "unique component",
            "present" if est.has_unique else "unavailable (NOT zero)")
        add(sec, "scope", "locus-level table context only — never drawn "
            "as positional signal")
        if result.estimate_report is not None:
            add(sec, "match", result.estimate_report.summary())

    if result.config.loci.subset_applied and result.config.loci.selected:
        add("subset", "applied subset",
            f"{len(result.config.loci.selected)} selected loci "
            f"(eligible population: {len(result.copies):,}); aggregates "
            "recomputed from original per-locus values, no reassignment")

    for i, w in enumerate(result.warnings, 1):
        add("warnings", f"warning {i}", w)

    return pd.DataFrame(rows, columns=["section", "metric", "value"])
