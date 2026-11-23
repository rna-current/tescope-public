"""Project/analysis configuration: the reproducibility contract.

A saved TEscope project file (YAML) fully specifies an analysis: assembly,
repeat, filters, tracks (with groups/replicates), binning, aggregation,
normalization, plotting and output options, plus the TEscope version that
wrote it. GUI and CLI both consume this model and call the same
:func:`tescope.core.analysis.run_analysis`.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field

from ..core.aggregation.aggregate import VALID_NORMALIZATIONS, VALID_STATISTICS
from ..core.filtering.filters import FilterConfig
from ..plotting.plotstyle import PlotStyleConfig


class ConfigError(ValueError):
    pass


@dataclass
class RepeatSelection:
    name: str  # RepeatMasker repeat name (subfamily), e.g. "L1HS"

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class TrackConfig:
    file: str
    name: str = ""  # defaults to basename
    kind: str = ""  # "continuous" | "interval" | "point"; "" = sniff
    format: str = ""  # "bigwig" | "bedgraph" | ... ; "" = sniff
    sample: str = ""
    group: str = ""
    replicate: int | None = None
    strand: str | None = None
    visible: bool = True
    color: str = ""  # display only
    # BAM-only options (validated: setting them on other inputs FAILS,
    # because they change the analyzed observation set)
    min_mapq: int = 0
    exclude_duplicates: bool = True
    mapping_class: str = ""  # "" / all | unique | multi (BAM/CRAM only)
    mapping_method: str = "nh"  # nh | mapq (aligner-dependent; see docs)
    unique_mapq: int = 30  # threshold for the mapq method

    def display_name(self) -> str:
        return self.name or os.path.basename(self.file)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class AnnotationConfig:
    """Consensus-feature annotation settings (metadata/visualization layer).

    Enabling annotation NEVER changes numerical analysis outputs; it
    attaches database (provider) and/or user-defined features to results
    and figures. ``user_features`` records are inline mappings with
    name/start/end (+ optional label/note), 0-based half-open on the
    consensus axis; ``user_features_file`` imports a TSV/YAML/JSON file.
    """

    enabled: bool = False
    provider: str = "dfam"  # "dfam" | "none"
    user_features: list = field(default_factory=list)  # inline records
    user_features_file: str = ""
    # imported feature files (docs/feature-import.md): list of mappings
    # {path, target: consensus|genome (optional for template TSV, which
    # declares it), reference: model name / assembly (optional when the
    # file declares it)}
    feature_files: list = field(default_factory=list)
    show_in_plots: bool = True
    required: bool = False  # fail the run if no annotation is available

    def validate(self) -> None:
        if self.provider not in ("dfam", "none", ""):
            raise ConfigError(f"annotations.provider must be 'dfam' or 'none', got {self.provider!r}")
        for i, rec in enumerate(self.user_features or []):
            if not isinstance(rec, dict):
                raise ConfigError(f"annotations.user_features[{i}] must be a mapping")
            for k in ("name", "start", "end"):
                if k not in rec:
                    raise ConfigError(f"annotations.user_features[{i}] is missing '{k}'")
        for i, rec in enumerate(self.feature_files or []):
            if not isinstance(rec, dict) or "path" not in rec:
                raise ConfigError(
                    f"annotations.feature_files[{i}] must be a mapping with 'path'"
                )
            unknown = set(rec) - {"path", "target", "reference"}
            if unknown:
                raise ConfigError(
                    f"annotations.feature_files[{i}]: unknown keys {sorted(unknown)}"
                )
            if rec.get("target", "") not in ("", "consensus", "genome"):
                raise ConfigError(
                    f"annotations.feature_files[{i}]: target must be "
                    "'consensus' or 'genome'"
                )

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class SequenceConfig:
    """Consensus-sequence source settings (display/metadata layer).

    Sequence display NEVER changes numerical analysis outputs, and
    sequence availability is never required to run an analysis: an
    unavailable/incompatible source degrades to an explained
    "unavailable" state. Sources: ``dfam`` (automatic, cached),
    ``library``/``fasta`` (user-provided local file — TEscope reads it in
    place and never copies or redistributes it), ``none``.
    """

    enabled: bool = False
    source: str = "dfam"  # "dfam" | "library" | "fasta" | "none"
    path: str = ""  # user-provided library/FASTA file (library/fasta sources)
    family: str = ""  # library record to use ("" = the analysis repeat name)
    library_name: str = ""  # optional user-entered provenance label
    show_track: bool = True  # draw the aligned sequence track when verified

    def validate(self) -> None:
        if self.source not in ("dfam", "library", "fasta", "none", ""):
            raise ConfigError(
                "sequence.source must be dfam/library/fasta/none, "
                f"got {self.source!r}"
            )
        if self.enabled and self.source in ("library", "fasta") and not self.path:
            raise ConfigError(
                f"sequence.source '{self.source}' requires sequence.path "
                "(the local FASTA/library file)"
            )

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class LociConfig:
    """Locus-level context: upstream estimates + reproducible subsets.

    ``estimates_file`` is a tescope-locus-estimates TSV (table-only
    context; never positional signal — docs/feature-import.md).
    ``selected`` is a list of stable locus IDs
    (``chrom:start-end(strand)``, 0-based half-open);
    ``subset_applied`` records whether the subset analysis was
    explicitly applied (browsing/selecting alone never changes any
    aggregate). Loading a pre-2 schema file simply means "no locus
    features in use" — absent fields are never reinterpreted.
    """

    estimates_file: str = ""
    selected: list[str] = field(default_factory=list)
    subset_applied: bool = False

    def validate(self) -> None:
        if self.subset_applied and len(self.selected) < 2:
            raise ConfigError(
                "loci.subset_applied requires >= 2 loci in loci.selected"
            )
        for s in self.selected:
            if not isinstance(s, str) or ":" not in s:
                raise ConfigError(
                    f"loci.selected entries must be stable locus IDs "
                    f"('chrom:start-end(strand)'), got {s!r}"
                )

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class AggregationConfig:
    bin_size: int = 10  # consensus bp per bin
    statistic: str = "mean"  # primary statistic for plots/group combining

    def validate(self) -> None:
        if self.bin_size < 1:
            raise ConfigError("aggregation.bin_size must be >= 1")
        if self.statistic not in VALID_STATISTICS:
            raise ConfigError(f"aggregation.statistic must be one of {VALID_STATISTICS}")


@dataclass
class PlotConfig:
    show_sem: bool = True
    y_axis: str = "independent"  # "independent" | "shared"
    normalization: str = "none"  # display-only: "none" | "zscore" | "minmax"
    show_replicates: bool = False  # plot individual replicates alongside group mean
    dpi: int = 300

    def validate(self) -> None:
        if self.y_axis not in ("independent", "shared"):
            raise ConfigError("plot.y_axis must be 'independent' or 'shared'")
        if self.normalization not in VALID_NORMALIZATIONS:
            raise ConfigError(f"plot.normalization must be one of {VALID_NORMALIZATIONS}")


@dataclass
class OutputConfig:
    dir: str = "tescope-results"
    table_formats: list[str] = field(default_factory=lambda: ["tsv"])
    figure_formats: list[str] = field(default_factory=lambda: ["png", "svg"])

    def validate(self) -> None:
        for f in self.table_formats:
            if f not in ("tsv", "csv"):
                raise ConfigError(f"output.table_formats: unknown format {f!r}")
        for f in self.figure_formats:
            if f not in ("png", "svg", "pdf", "tiff"):
                raise ConfigError(f"output.figure_formats: unknown format {f!r}")


@dataclass
class AnalysisConfig:
    assembly: str
    repeat: RepeatSelection
    tracks: list[TrackConfig] = field(default_factory=list)
    filters: FilterConfig = field(default_factory=FilterConfig)
    aggregation: AggregationConfig = field(default_factory=AggregationConfig)
    plot: PlotConfig = field(default_factory=PlotConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    annotations: AnnotationConfig = field(default_factory=AnnotationConfig)
    sequence: SequenceConfig = field(default_factory=SequenceConfig)
    loci: LociConfig = field(default_factory=LociConfig)
    # pure rendering layer — never read by analysis code (see plotstyle.py)
    style: PlotStyleConfig = field(default_factory=PlotStyleConfig)
    chrom_alias: bool = True  # allow unambiguous chr1<->1 / chrM<->MT reconciliation
    tescope_version: str = ""
    annotation_provenance: dict = field(default_factory=dict)  # filled at run time

    def validate(self) -> None:
        if not self.assembly:
            raise ConfigError("assembly is required")
        if not self.repeat.name:
            raise ConfigError("repeat.name is required")
        if not self.tracks:
            raise ConfigError("at least one track is required")
        seen = set()
        for t in self.tracks:
            if t.kind and t.kind not in ("continuous", "interval", "point"):
                raise ConfigError(
                    f"track {t.display_name()}: kind must be continuous/interval/point"
                )
            if t.mapping_class not in ("", "all", "unique", "multi"):
                raise ConfigError(
                    f"track {t.display_name()}: mapping_class must be "
                    f"all/unique/multi, got {t.mapping_class!r}"
                )
            if t.mapping_method not in ("nh", "mapq"):
                raise ConfigError(
                    f"track {t.display_name()}: mapping_method must be 'nh' or 'mapq'"
                )
            if t.mapping_class in ("unique", "multi"):
                is_bam = (t.format == "bam"
                          or t.file.lower().endswith((".bam", ".cram")))
                if not is_bam:
                    raise ConfigError(
                        f"Requested mapping view unavailable: track "
                        f"'{t.display_name()}' is pre-aggregated and contains "
                        "no read-level mapping information (mapping "
                        "classification requires BAM/CRAM)."
                    )
                if t.mapping_class == "multi" and t.mapping_method == "mapq":
                    raise ConfigError(
                        f"track {t.display_name()}: mapping_class 'multi' "
                        "requires mapping_method 'nh' (MAPQ cannot positively "
                        "identify multi-mappers)."
                    )
            key = (t.display_name(), t.kind)
            if key in seen:
                raise ConfigError(f"duplicate track name {t.display_name()!r}")
            seen.add(key)
        self.aggregation.validate()
        self.plot.validate()
        self.output.validate()
        self.annotations.validate()
        self.sequence.validate()
        self.loci.validate()
        try:
            self.style.validate()
        except ValueError as e:
            raise ConfigError(str(e)) from e

    # -- (de)serialization -------------------------------------------------

    def to_dict(self) -> dict:
        from tescope import __version__

        return {
            "tescope_version": self.tescope_version or __version__,
            # explicit schema version (1 = rc2-era files without it);
            # absent new sections in old files mean "feature not used",
            # never a reinterpretation
            "config_schema": 2,
            "assembly": self.assembly,
            "repeat": self.repeat.to_dict(),
            "filters": self.filters.to_dict(),
            "tracks": [t.to_dict() for t in self.tracks],
            "aggregation": asdict(self.aggregation),
            "plot": asdict(self.plot),
            "output": asdict(self.output),
            "annotations": self.annotations.to_dict(),
            "sequence": self.sequence.to_dict(),
            "loci": self.loci.to_dict(),
            "style": self.style.to_dict(),
            "chrom_alias": self.chrom_alias,
            "annotation_provenance": self.annotation_provenance,
        }

    @classmethod
    def from_dict(cls, d: dict) -> AnalysisConfig:
        try:
            repeat_d = d["repeat"]
            repeat = RepeatSelection(name=str(repeat_d["name"]))
        except (KeyError, TypeError):
            raise ConfigError("config must contain repeat: {name: <repeat name>}") from None
        if "assembly" not in d:
            raise ConfigError("config must contain 'assembly'")
        tracks = []
        for i, td in enumerate(d.get("tracks") or []):
            if "file" not in td:
                raise ConfigError(f"tracks[{i}] is missing 'file'")
            known = {k: v for k, v in td.items() if k in TrackConfig.__dataclass_fields__}
            unknown = set(td) - set(known)
            if unknown:
                raise ConfigError(f"tracks[{i}] has unknown keys: {sorted(unknown)}")
            tracks.append(TrackConfig(**known))
        cfg = cls(
            assembly=str(d["assembly"]),
            repeat=repeat,
            tracks=tracks,
            filters=FilterConfig.from_dict(d.get("filters") or {}),
            aggregation=_sub(AggregationConfig, d.get("aggregation")),
            plot=_sub(PlotConfig, d.get("plot")),
            output=_sub(OutputConfig, d.get("output")),
            annotations=_sub(AnnotationConfig, d.get("annotations")),
            sequence=_sub(SequenceConfig, d.get("sequence")),
            loci=_sub(LociConfig, d.get("loci")),
            style=PlotStyleConfig.from_dict(d.get("style") or {}),
            chrom_alias=bool(d.get("chrom_alias", True)),
            tescope_version=str(d.get("tescope_version", "")),
            annotation_provenance=dict(d.get("annotation_provenance") or {}),
        )
        cfg.validate()
        return cfg


def _sub(cls, d: dict | None):
    d = d or {}
    known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
    unknown = set(d) - set(known)
    if unknown:
        raise ConfigError(f"{cls.__name__.lower()}: unknown keys {sorted(unknown)}")
    return cls(**known)
