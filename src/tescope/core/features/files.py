"""Feature-file import: TE features on consensus OR genomic coordinates.

Formats (docs/feature-import.md):

- ``tescope-features v1`` TSV: header directives declare target
  (consensus|genome), reference, and the coordinate convention
  (always 0-based half-open for the template).
- BED (.bed[.gz]): genome target, 0-based half-open **by the BED
  specification** — never guessed from values.
- GFF3 (.gff3/.gff[.gz]): genome target, 1-based inclusive **by the
  GFF3 specification**, converted internally; ``ID``/``Parent``
  attributes are used for identity/hierarchy.

Internal convention: 0-based half-open everywhere; the original
convention is recorded. Amino-acid coordinates are rejected with an
explanation (no x3 conversion, no frame guessing). Validation is
collected into a report — a file either imports fully validated or
fails with every problem listed.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote

from .model import SOURCE_USER, ConsensusFeature

_MAGIC = "tescope-features"
_HEADER_RE = re.compile(r"^#\s*([A-Za-z_]+)\s*:\s*(.*?)\s*$")
_TSV_COLUMNS = (
    "feature_id", "start", "end", "strand", "label", "type",
    "parent_id", "source",
)
_VALID_STRANDS = ("+", "-", ".")


class FeatureFileError(ValueError):
    """User-facing import failure; message lists every validation error."""


@dataclass
class RichFeature:
    """One imported feature, internally 0-based half-open."""

    feature_id: str
    start: int
    end: int
    strand: str = "."
    label: str = ""
    type: str = ""
    parent_id: str = ""
    source: str = ""
    chrom: str = ""  # genome-target features only

    def length(self) -> int:
        return self.end - self.start


@dataclass
class FeatureFileResult:
    path: str
    target: str  # "consensus" | "genome"
    reference: str  # consensus model name or assembly, as declared
    source_convention: str  # recorded verbatim
    format: str  # "tsv" | "bed" | "gff3"
    features: list[RichFeature] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def consensus_features(self) -> list[ConsensusFeature]:
        """Convert to the existing user-feature lane (consensus target)."""
        if self.target != "consensus":
            raise FeatureFileError(
                f"{self.path}: genome-target features cannot be attached to "
                "the consensus axis"
            )
        out = []
        for f in self.features:
            note = f.source or f"imported from {Path(self.path).name}"
            if f.parent_id:
                note += f"; parent={f.parent_id}"
            out.append(
                ConsensusFeature(
                    name=f.label or f.feature_id,
                    start=f.start,
                    end=f.end,
                    origin=SOURCE_USER,
                    type=f.type,
                    note=note,
                )
            )
        return out


@dataclass
class GenomicFeatureSet:
    """Genome-target features attached to an analysis (one assembly)."""

    assembly: str
    path: str
    format: str
    features: list[RichFeature] = field(default_factory=list)
    status: str = ""

    def overlapping(
        self, chrom: str, start: int, end: int, flank: int = 0
    ) -> list[RichFeature]:
        a, b = start - flank, end + flank
        return [
            f
            for f in self.features
            if f.chrom == chrom and f.start < b and f.end > a
        ]


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def parse_feature_file(path: str) -> FeatureFileResult:
    p = Path(path)
    if not p.exists():
        raise FeatureFileError(f"feature file not found: {p}")
    name = p.name.lower()
    base = name[:-3] if name.endswith(".gz") else name
    if base.endswith((".gff3", ".gff")):
        return _parse_gff3(p)
    if base.endswith(".bed"):
        return _parse_bed(p)
    return _parse_template_tsv(p)


def _finish(result: FeatureFileResult, errors: list[str]) -> FeatureFileResult:
    errors.extend(_validate(result))
    if errors:
        shown = "\n  ".join(errors[:12])
        more = f"\n  (+{len(errors) - 12} more)" if len(errors) > 12 else ""
        raise FeatureFileError(
            f"{result.path}: {len(errors)} validation error(s):\n  {shown}{more}"
        )
    return result


def _validate(result: FeatureFileResult) -> list[str]:
    errors: list[str] = []
    seen: set[str] = set()
    ids = {f.feature_id for f in result.features}
    for f in result.features:
        where = f"feature '{f.feature_id}'"
        if not f.feature_id:
            errors.append("feature with empty feature_id")
            continue
        if f.feature_id in seen:
            errors.append(f"duplicate feature_id '{f.feature_id}'")
        seen.add(f.feature_id)
        if f.end <= f.start:
            errors.append(f"{where}: end ({f.end}) must be > start ({f.start})")
        if f.start < 0:
            errors.append(f"{where}: negative start")
        if f.strand not in _VALID_STRANDS:
            errors.append(f"{where}: strand must be one of {_VALID_STRANDS}")
        if f.parent_id and f.parent_id not in ids:
            errors.append(
                f"{where}: parent_id '{f.parent_id}' does not match any "
                "feature_id in this file"
            )
        if result.target == "genome" and not f.chrom:
            errors.append(f"{where}: genome-target feature without a contig")
    if not result.features:
        errors.append("no features found")
    return errors


def _parse_template_tsv(p: Path) -> FeatureFileResult:
    with _open_text(p) as fh:
        first = fh.readline()
        if _MAGIC not in first:
            raise FeatureFileError(
                f"{p}: unrecognized format — expected a '{_MAGIC} v1' TSV, "
                "a .bed, or a .gff3 file (see docs/feature-import.md)"
            )
        directives: dict[str, str] = {}
        header: list[str] | None = None
        rows: list[tuple[int, list[str]]] = []
        for lineno, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            m = _HEADER_RE.match(line)
            if m and header is None:
                directives[m.group(1).lower()] = m.group(2)
                continue
            if line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                continue
            rows.append((lineno, line.split("\t")))

    units = directives.get("units", "").strip().lower()
    if units in ("aa", "amino-acids", "aminoacids", "codons"):
        raise FeatureFileError(
            f"{p}: amino-acid coordinates are not supported. TEscope does "
            "not multiply by three or guess CDS frames; convert to "
            "nucleotide coordinates yourself if you know the frame."
        )
    target = directives.get("target", "").strip().lower()
    if target not in ("consensus", "genome"):
        raise FeatureFileError(
            f"{p}: '# target:' directive must be 'consensus' or 'genome'"
        )
    convention = directives.get("coordinates", "0-based-half-open").strip()
    if convention != "0-based-half-open":
        raise FeatureFileError(
            f"{p}: the TSV template is defined as 0-based-half-open; "
            f"got {convention!r}. Use BED/GFF3 for their native conventions."
        )
    if header is None:
        raise FeatureFileError(f"{p}: missing column header line")
    unknown = set(h.strip() for h in header) - set(_TSV_COLUMNS)
    if unknown:
        raise FeatureFileError(f"{p}: unknown column(s) {sorted(unknown)}")
    for req in ("feature_id", "start", "end"):
        if req not in header:
            raise FeatureFileError(f"{p}: missing required column '{req}'")

    result = FeatureFileResult(
        path=str(p),
        target=target,
        reference=directives.get("reference", ""),
        source_convention="0-based half-open (template)",
        format="tsv",
    )
    errors: list[str] = []
    idx = {h: i for i, h in enumerate(header)}

    def get(vals: list[str], col: str) -> str:
        i = idx.get(col)
        return vals[i].strip() if i is not None and i < len(vals) else ""

    for lineno, vals in rows:
        try:
            start = int(get(vals, "start"))
            end = int(get(vals, "end"))
        except ValueError:
            errors.append(f"line {lineno}: start/end must be integers")
            continue
        chrom = ""
        fid = get(vals, "feature_id")
        if target == "genome":
            # genomic template IDs may carry the contig as 'chrom:fid';
            # otherwise a 'source'-declared contig column is not part of
            # the v1 template — require the chrom:fid form.
            if ":" in fid:
                chrom, fid = fid.split(":", 1)
            else:
                errors.append(
                    f"line {lineno}: genome-target feature_id must be "
                    "'contig:feature_id' in the v1 template (or use BED/GFF3)"
                )
        result.features.append(
            RichFeature(
                feature_id=fid,
                start=start,
                end=end,
                strand=get(vals, "strand") or ".",
                label=get(vals, "label"),
                type=get(vals, "type"),
                parent_id=get(vals, "parent_id"),
                source=get(vals, "source"),
                chrom=chrom,
            )
        )
    return _finish(result, errors)


def _parse_bed(p: Path) -> FeatureFileResult:
    result = FeatureFileResult(
        path=str(p),
        target="genome",
        reference="",  # assembly comes from the config entry
        source_convention="0-based half-open (BED specification)",
        format="bed",
    )
    errors: list[str] = []
    with _open_text(p) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 3:
                errors.append(f"line {lineno}: BED needs >= 3 fields")
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                errors.append(f"line {lineno}: non-integer coordinates")
                continue
            name = f[3] if len(f) > 3 and f[3] not in (".", "") else f"{p.stem}:{lineno}"
            strand = f[5] if len(f) > 5 and f[5] in ("+", "-") else "."
            result.features.append(
                RichFeature(
                    feature_id=name,
                    start=start,
                    end=end,
                    strand=strand,
                    label=name,
                    type="BED feature",
                    chrom=f[0],
                )
            )
    # BED names may repeat; qualify duplicates instead of failing
    seen: dict[str, int] = {}
    for f in result.features:
        n = seen.get(f.feature_id, 0)
        seen[f.feature_id] = n + 1
        if n:
            f.feature_id = f"{f.feature_id}#{n + 1}"
    if any(v > 1 for v in seen.values()):
        result.warnings.append(
            "duplicate BED names were suffixed with #2, #3, … to keep "
            "feature identities distinct (records were NOT merged)"
        )
    return _finish(result, errors)


def _parse_gff3(p: Path) -> FeatureFileResult:
    result = FeatureFileResult(
        path=str(p),
        target="genome",
        reference="",
        source_convention="1-based inclusive (GFF3 specification)",
        format="gff3",
    )
    errors: list[str] = []
    auto = 0
    with _open_text(p) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) != 9:
                errors.append(f"line {lineno}: GFF3 needs 9 fields")
                continue
            try:
                start1, end1 = int(f[3]), int(f[4])
            except ValueError:
                errors.append(f"line {lineno}: non-integer coordinates")
                continue
            attrs = {}
            for kv in f[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k.strip()] = unquote(v.strip())
            fid = attrs.get("ID", "")
            if not fid:
                auto += 1
                fid = f"{p.stem}:auto{auto}"
            result.features.append(
                RichFeature(
                    feature_id=fid,
                    start=start1 - 1,  # GFF3 1-based inclusive -> internal
                    end=end1,
                    strand=f[6] if f[6] in ("+", "-") else ".",
                    label=attrs.get("Name", fid),
                    type=f[2],
                    parent_id=attrs.get("Parent", ""),
                    source=f[1] if f[1] != "." else "",
                    chrom=f[0],
                )
            )
    return _finish(result, errors)


def attach_genome_features(
    parsed: FeatureFileResult,
    assembly: str,
    declared_assembly: str,
    known_chroms: list[str] | None,
) -> GenomicFeatureSet:
    """Validate a genome-target file against the analysis assembly."""
    if parsed.target != "genome":
        raise FeatureFileError(f"{parsed.path}: not a genome-target file")
    ref = parsed.reference or declared_assembly
    if not ref:
        raise FeatureFileError(
            f"{parsed.path}: no assembly declared — set 'reference' in the "
            "file header or 'reference' on the config entry"
        )
    if ref != assembly:
        raise FeatureFileError(
            f"{parsed.path}: declares assembly {ref!r} but the analysis "
            f"uses {assembly!r}; explicit conversion is required "
            "(no silent aliasing)"
        )
    if known_chroms:
        known = set(known_chroms)
        bad = sorted({f.chrom for f in parsed.features} - known)
        if bad:
            raise FeatureFileError(
                f"{parsed.path}: contig(s) not in {assembly}: "
                f"{', '.join(bad[:5])}"
                + (f" (+{len(bad) - 5} more)" if len(bad) > 5 else "")
            )
    return GenomicFeatureSet(
        assembly=assembly,
        path=parsed.path,
        format=parsed.format,
        features=parsed.features,
        status=(
            f"{len(parsed.features)} genomic feature(s) from "
            f"{Path(parsed.path).name} ({parsed.source_convention})"
        ),
    )
