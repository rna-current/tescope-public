"""Analysis orchestrator: config in, results out.

This is the single entry point shared by the CLI and the GUI. It performs
no plotting and no terminal/window I/O; progress and log messages go
through optional callbacks.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..config.model import AnalysisConfig, TrackConfig
from .aggregation.aggregate import aggregate_profile, combine_replicates
from .annotation.schema import add_derived_columns, drop_inconsistent_hits
from .chromnames import ChromReport, reconcile
from .coordinates.mapping import make_bins
from .features.model import FeatureSet
from .features.resolve import AnnotationRequiredError, resolve_features
from .filtering.filters import FilterReport, apply_filters, classify_completeness
from .projection.continuous import CopyBinMatrix, project_continuous
from .projection.coverage import element_coverage
from .projection.intervals import interval_occurrence, project_intervals, project_points
from .reference.store import ReferenceStore
from .sequence.model import SequenceResolution
from .sequence.resolve import resolve_sequence
from .tracks.base import ContinuousTrack, IntervalTrack, PointTrack, TrackKind
from .tracks.loader import open_track

ProgressCallback = Callable[[str, int, int], None]  # (message, step, total_steps)


class AnalysisError(RuntimeError):
    """User-facing analysis failure with an actionable message."""


@dataclass
class TrackResult:
    config: TrackConfig
    kind: TrackKind
    chrom_report: ChromReport | None = None
    # continuous
    matrix: CopyBinMatrix | None = None
    profile: pd.DataFrame | None = None
    # interval
    projected: pd.DataFrame | None = None
    occurrence: pd.DataFrame | None = None
    # point
    points: pd.DataFrame | None = None
    point_occurrence: pd.DataFrame | None = None
    n_source_items: int | None = None

    @property
    def name(self) -> str:
        return self.config.display_name()


@dataclass
class AnalysisResult:
    config: AnalysisConfig
    repeat_name: str
    cons_length: int
    bin_edges: np.ndarray
    n_annotated: int
    filter_report: FilterReport
    copies: pd.DataFrame
    element_covered: np.ndarray  # geometry-only (n_copies, n_bins)
    tracks: list[TrackResult] = field(default_factory=list)
    group_profiles: dict[str, pd.DataFrame] = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    # consensus-feature annotation: metadata/visualization layer ONLY —
    # never participates in projection/aggregation/normalization
    features: FeatureSet | None = None
    # consensus sequence: display/metadata layer ONLY, same contract;
    # attached sequence is guaranteed model-compatible (verified)
    sequence: SequenceResolution | None = None
    # genome-target user feature sets (context layer; per-locus display)
    genome_features: list = field(default_factory=list)
    # upstream per-locus estimates: table-only context (never positional)
    locus_estimates: object | None = None
    estimate_report: object | None = None
    # lazily built per-locus table (core.loci.table.build_locus_table)
    _locus_table: pd.DataFrame | None = None

    @property
    def n_selected(self) -> int:
        return len(self.copies)

    def locus_table(self) -> pd.DataFrame:
        """Per-locus table (cached). Derived from existing data only."""
        if self._locus_table is None:
            from .loci.estimates import match_estimates
            from .loci.table import build_locus_table

            table = build_locus_table(self)
            if self.locus_estimates is not None:
                table, self.estimate_report = match_estimates(
                    self.locus_estimates, table, self.config.assembly
                )
            self._locus_table = table
        return self._locus_table

    def continuous_tracks(self) -> list[TrackResult]:
        return [t for t in self.tracks if t.kind is TrackKind.CONTINUOUS]

    def interval_tracks(self) -> list[TrackResult]:
        return [t for t in self.tracks if t.kind is TrackKind.INTERVAL]

    def point_tracks(self) -> list[TrackResult]:
        return [t for t in self.tracks if t.kind is TrackKind.POINT]


def run_analysis(
    config: AnalysisConfig,
    store: ReferenceStore | None = None,
    progress: ProgressCallback | None = None,
) -> AnalysisResult:
    config.validate()
    store = store or ReferenceStore()
    total_steps = 3 + len(config.tracks)
    step = 0

    def tick(msg: str):
        nonlocal step
        step += 1
        if progress:
            progress(msg, step, total_steps)

    tick(f"Loading {config.assembly} annotation")
    repeat_name = config.repeat.name
    ann = store.select_repeat(config.assembly, repeat_name)
    cons_lengths = store.consensus_lengths(config.assembly)
    try:
        cons_length = int(cons_lengths[repeat_name])
    except (KeyError, TypeError):
        cons_length = 0
    if cons_length <= 0:
        raise AnalysisError(
            f"Consensus model length for '{repeat_name}' could not be determined from "
            "the annotation (missing repeat-consensus coordinates). TEscope does not "
            "silently fall back to normalized interval mode."
        )
    prov = store.provenance(config.assembly).to_dict()
    config.annotation_provenance = prov

    tick("Filtering TE copies")
    warnings: list[str] = []
    ann, n_inc = drop_inconsistent_hits(ann, cons_length)
    if n_inc:
        warnings.append(
            f"{n_inc} {repeat_name} hits excluded: consensus coordinates exceed "
            f"the model length ({cons_length} bp) — annotation artifact."
        )
        if len(ann) == 0:
            raise AnalysisError(
                f"All {repeat_name} hits have consensus coordinates inconsistent "
                "with the derived model length; the annotation may be corrupt."
            )
    ann = add_derived_columns(ann, cons_lengths)
    ann["completeness"] = classify_completeness(
        ann,
        end_tolerance_frac=config.filters.end_tolerance_frac,
        end_tolerance_bp=config.filters.end_tolerance_bp,
    )
    copies, filter_report = apply_filters(ann, config.filters)
    if len(copies) == 0:
        raise AnalysisError(
            "No TE annotations satisfy the current filters. "
            f"({filter_report.n_input} {repeat_name} copies before filtering.)"
        )

    if len(copies) < 10:
        warnings.append(
            f"Only {len(copies)} TE copies remain after filtering; "
            "aggregate profiles will be noisy."
        )

    bin_edges = make_bins(cons_length, config.aggregation.bin_size)
    elem_cov = element_coverage(copies, bin_edges)

    ref_chroms = _reference_chroms(store, config.assembly, ann)

    result = AnalysisResult(
        config=config,
        repeat_name=repeat_name,
        cons_length=cons_length,
        bin_edges=bin_edges,
        n_annotated=filter_report.n_input,
        filter_report=filter_report,
        copies=copies,
        element_covered=elem_cov,
        provenance=prov,
        warnings=warnings,
    )

    for tcfg in config.tracks:
        tick(f"Projecting {tcfg.display_name()}")
        tres = _run_track(tcfg, copies, cons_length, bin_edges, elem_cov, ref_chroms, config)
        result.tracks.append(tres)
        if tres.chrom_report and tres.chrom_report.aliased:
            warnings.append(
                f"Track '{tres.name}': chromosome names reconciled via alias "
                f"({tres.chrom_report.summary()})."
            )
        if tres.matrix and tres.matrix.n_skipped_chrom:
            warnings.append(
                f"Track '{tres.name}': {tres.matrix.n_skipped_chrom} TE copies skipped "
                "(chromosome absent from track)."
            )

    tick("Combining groups")
    _combine_groups(result, config)

    if config.annotations.enabled:
        # strictly after all numerical work; failures degrade to a warning
        try:
            result.features = resolve_features(
                config.annotations, repeat_name, cons_length,
                cache_dir=store.cache_dir,
            )
        except AnnotationRequiredError:
            raise
        except Exception as e:  # noqa: BLE001 - annotation must not fail analysis
            warnings.append(f"Annotation unavailable: {e}")
        else:
            if result.features.status:
                warnings.append(f"Annotation: {result.features.status}")

    # imported feature files (explicit user inputs: malformed or
    # incompatible files FAIL the run with the full validation report —
    # silent degradation would hide a mistake the user must fix)
    if config.annotations.feature_files:
        _attach_feature_files(result, config, store)

    # upstream per-locus estimates: table-only context layer
    if config.loci.estimates_file:
        from .loci.estimates import EstimateError, load_estimates

        try:
            result.locus_estimates = load_estimates(config.loci.estimates_file)
        except (OSError, EstimateError) as e:
            raise AnalysisError(f"Locus-estimate import failed: {e}") from e
        est = result.locus_estimates
        warnings.append(
            f"Locus estimates attached ({est.assay}, units: {est.units}"
            + (f", tool: {est.tool}" if est.tool else "")
            + ") — locus-level table context only, never positional signal."
        )

    if config.sequence.enabled:
        # display/metadata only — resolution failure NEVER fails analysis
        try:
            result.sequence = resolve_sequence(
                config.sequence, repeat_name, cons_length,
                cache_dir=store.cache_dir,
            )
        except Exception as e:  # noqa: BLE001 - sequence must not fail analysis
            result.sequence = SequenceResolution(
                detail=f"sequence resolution failed: {e}"
            )
        if not result.sequence.verified:
            warnings.append(
                f"Consensus sequence {result.sequence.status}: "
                f"{result.sequence.detail}"
            )
    return result


def _reference_chroms(store: ReferenceStore, assembly: str, ann: pd.DataFrame) -> list[str]:
    sizes = store.chrom_sizes(assembly)
    if sizes:
        return list(sizes)
    return [str(c) for c in ann["chrom"].unique()]


def _open_configured_track(tcfg: TrackConfig):
    kwargs: dict = {}
    if tcfg.format == "bam" or tcfg.file.lower().endswith((".bam", ".cram")):
        kwargs.update(
            min_mapq=tcfg.min_mapq,
            exclude_duplicates=tcfg.exclude_duplicates,
            strand=tcfg.strand,
            mapping_class=tcfg.mapping_class,
            mapping_method=tcfg.mapping_method,
            unique_mapq=tcfg.unique_mapq,
        )
    else:
        kwargs.update(
            sample=tcfg.sample,
            group=tcfg.group,
            replicate=tcfg.replicate,
        )
    try:
        return open_track(
            tcfg.file,
            name=tcfg.name or None,
            kind=tcfg.kind or None,
            format=tcfg.format or None,
            **kwargs,
        )
    except FileNotFoundError:
        raise AnalysisError(f"Track file not found: {tcfg.file}") from None
    except (ValueError, RuntimeError) as e:
        raise AnalysisError(f"Track '{tcfg.display_name()}': {e}") from e


def _run_track(
    tcfg: TrackConfig,
    copies: pd.DataFrame,
    cons_length: int,
    bin_edges: np.ndarray,
    elem_cov: np.ndarray,
    ref_chroms: list[str],
    config: AnalysisConfig,
) -> TrackResult:
    track = _open_configured_track(tcfg)
    try:
        report = reconcile(
            list(track.chroms()), ref_chroms, allow_alias=config.chrom_alias
        )
        copy_chroms = set(copies["chrom"].astype(str))
        translatable = {c for c in copy_chroms if report.translate(c)}
        if not translatable:
            raise AnalysisError(
                f"Track '{tcfg.display_name()}' shares no chromosomes with the "
                f"selected TE copies on {config.assembly} "
                f"(track chromosomes e.g. {', '.join(list(track.chroms())[:3])}). "
                "Check that the track was produced on the selected assembly."
            )
        if isinstance(track, ContinuousTrack):
            matrix = project_continuous(
                track, copies, cons_length, config.aggregation.bin_size, report
            )
            profile = aggregate_profile(matrix)
            return TrackResult(
                config=tcfg,
                kind=TrackKind.CONTINUOUS,
                chrom_report=report,
                matrix=matrix,
                profile=profile,
            )
        if isinstance(track, IntervalTrack):
            projected = project_intervals(track, copies, report)
            occurrence = interval_occurrence(
                projected, elem_cov, copies.index.to_numpy(), bin_edges
            )
            return TrackResult(
                config=tcfg,
                kind=TrackKind.INTERVAL,
                chrom_report=report,
                projected=projected,
                occurrence=occurrence,
                n_source_items=track.count(),
            )
        if isinstance(track, PointTrack):
            points = project_points(track, copies, report)
            point_occ = _point_occurrence(points, elem_cov, bin_edges)
            return TrackResult(
                config=tcfg,
                kind=TrackKind.POINT,
                chrom_report=report,
                points=points,
                point_occurrence=point_occ,
                n_source_items=track.count(),
            )
        raise AnalysisError(f"Track '{tcfg.display_name()}': unsupported track class")
    finally:
        track.close()


def _attach_feature_files(
    result: AnalysisResult, config: AnalysisConfig, store: ReferenceStore
) -> None:
    from .features.files import (
        FeatureFileError,
        attach_genome_features,
        parse_feature_file,
    )
    from .features.model import FeatureSet

    known_chroms = list(store.chrom_sizes(config.assembly) or [])
    for entry in config.annotations.feature_files:
        path = entry["path"]
        try:
            parsed = parse_feature_file(path)
            declared = entry.get("target", "")
            if declared and parsed.target != declared:
                raise FeatureFileError(
                    f"{path}: config declares target {declared!r} but the "
                    f"file is {parsed.target!r}-target"
                )
            if parsed.target == "consensus":
                ref = parsed.reference or entry.get("reference", "")
                if ref and ref != result.repeat_name:
                    raise FeatureFileError(
                        f"{path}: declares consensus model {ref!r} but the "
                        f"analysis model is {result.repeat_name!r} — "
                        "features withheld (name compatibility is required; "
                        "it is still not proof of sequence identity)"
                    )
                feats = parsed.consensus_features()
                for f in feats:
                    f.validate(result.cons_length)
                if result.features is None:
                    result.features = FeatureSet(
                        repeat_name=result.repeat_name,
                        cons_length=result.cons_length,
                    )
                result.features.user_features.extend(feats)
                result.warnings.append(
                    f"{len(feats)} consensus feature(s) imported from "
                    f"{path} (user-defined lane)."
                )
            else:
                gset = attach_genome_features(
                    parsed,
                    config.assembly,
                    entry.get("reference", ""),
                    known_chroms,
                )
                result.genome_features.append(gset)
                result.warnings.append(gset.status)
            for w in parsed.warnings:
                result.warnings.append(f"{path}: {w}")
        except (FeatureFileError, ValueError) as e:
            raise AnalysisError(f"Feature import failed: {e}") from e
        except OSError as e:
            raise AnalysisError(f"Feature import failed: {path}: {e}") from e


def _point_occurrence(
    points: pd.DataFrame, elem_cov: np.ndarray, bin_edges: np.ndarray
) -> pd.DataFrame:
    counts, _ = np.histogram(points["cons_position"].to_numpy(), bins=bin_edges)
    n_elements = elem_cov.sum(axis=0).astype("int64")
    with np.errstate(invalid="ignore", divide="ignore"):
        per_element = np.where(n_elements > 0, counts / n_elements, np.nan)
    return pd.DataFrame(
        {
            "consensus_bin": np.arange(len(bin_edges) - 1),
            "consensus_start": bin_edges[:-1],
            "consensus_end": bin_edges[1:],
            "point_count": counts.astype("int64"),
            "n_elements": n_elements,
            "points_per_element": per_element,
        }
    )


def _combine_groups(result: AnalysisResult, config: AnalysisConfig) -> None:
    groups: dict[str, list[pd.DataFrame]] = {}
    for t in result.continuous_tracks():
        if t.config.group and t.profile is not None:
            groups.setdefault(t.config.group, []).append(t.profile)
    for gname, profiles in groups.items():
        if len(profiles) >= 1:
            result.group_profiles[gname] = combine_replicates(
                profiles, statistic=config.aggregation.statistic
            )
