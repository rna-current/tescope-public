"""BAM/CRAM coverage as a convenience continuous track (optional, pysam).

Scope guard: TEscope derives simple per-base coverage from an indexed,
coordinate-sorted BAM/CRAM. It does **not** attempt multi-mapping read
rescue or TE expression quantification; reads in young TE families are
frequently multi-mapping and coverage there under-represents true signal.
This is stated in the docs and repeated in the QC output whenever a BAM
track is used.

Mapping classification (explicit, aligner-dependent — documented):
- method "nh": NH tag == 1 -> unique; NH > 1 -> multi; missing NH ->
  ambiguous. Requires the aligner to have written NH tags.
- method "mapq": MAPQ >= unique_mapq -> unique; below -> ambiguous.
  MAPQ alone cannot positively identify multi-mappers, so a "multi"
  selection is NOT available under the mapq method.
``mapping_class`` selects which classified reads contribute coverage:
"" / "all" (default, previous behavior), "unique", or "multi" (nh only).
Secondary/supplementary alignments remain excluded, so a multi-mapping
read is counted once (at its primary alignment) — never multiplied.
These options change the observation set: they are ANALYSIS settings,
recorded in config/provenance/QC, never silent.
"""

from __future__ import annotations

import os

import numpy as np

from .base import ContinuousTrack, TrackKind, TrackMeta

try:  # pragma: no cover - import probing
    import pysam  # type: ignore
except ImportError:  # pragma: no cover
    pysam = None


class BamCoverageError(RuntimeError):
    pass


def bam_available() -> bool:
    return pysam is not None


class BamCoverageTrack(ContinuousTrack):
    def __init__(
        self,
        meta: TrackMeta,
        *,
        min_mapq: int = 0,
        exclude_duplicates: bool = True,
        strand: str | None = None,  # None = both, '+' or '-' = only that strand
        create_index: bool = False,
        mapping_class: str = "",  # "" / "all" | "unique" | "multi"
        mapping_method: str = "nh",  # "nh" | "mapq"
        unique_mapq: int = 30,
    ):
        super().__init__(meta)
        if mapping_class not in ("", "all", "unique", "multi"):
            raise BamCoverageError(
                f"mapping_class must be all/unique/multi, got {mapping_class!r}")
        if mapping_method not in ("nh", "mapq"):
            raise BamCoverageError(
                f"mapping_method must be 'nh' or 'mapq', got {mapping_method!r}")
        if mapping_class == "multi" and mapping_method == "mapq":
            raise BamCoverageError(
                "mapping_class 'multi' requires the NH method: MAPQ alone "
                "cannot positively identify multi-mapping reads")
        self.mapping_class = mapping_class or "all"
        self.mapping_method = mapping_method
        self.unique_mapq = int(unique_mapq)
        if pysam is None:
            raise BamCoverageError(
                "BAM/CRAM support requires pysam (install tescope[bam]); "
                "pysam publishes no Windows wheels — use bigWig coverage on Windows"
            )
        path = meta.path
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        self.min_mapq = int(min_mapq)
        self.exclude_duplicates = bool(exclude_duplicates)
        self.strand = strand
        mode = "rc" if path.lower().endswith(".cram") else "rb"
        self._af = pysam.AlignmentFile(path, mode)
        if not self._af.has_index():
            if create_index:
                self._af.close()
                pysam.index(path)
                self._af = pysam.AlignmentFile(path, mode)
            else:
                raise BamCoverageError(
                    f"{path}: BAM index missing. Sort and index the file "
                    "(samtools sort + samtools index) or enable automatic indexing."
                )
        if self.mapping_class != "all" and self.mapping_method == "nh":
            # verify the actual file carries NH tags — fail clearly rather
            # than silently analyzing an unclassifiable observation set
            from .capabilities import _bam_probe

            _, has_nh = _bam_probe(path)
            if not has_nh:
                raise BamCoverageError(
                    f"{path}: mapping_class '{self.mapping_class}' with the NH "
                    "method requires NH tags, but this alignment file does not "
                    "contain them. Use mapping_method: mapq (unique only) or "
                    "re-align with an NH-writing aligner."
                )
        meta.metadata.setdefault("min_mapq", self.min_mapq)
        meta.metadata.setdefault("exclude_duplicates", self.exclude_duplicates)
        meta.metadata.setdefault("mapping_class", self.mapping_class)
        meta.metadata.setdefault("mapping_method", self.mapping_method)
        if self.mapping_method == "mapq":
            meta.metadata.setdefault("unique_mapq", self.unique_mapq)
        meta.metadata.setdefault(
            "caveat",
            "BAM-derived coverage; multi-mapping reads in TE regions are not rescued",
        )

    def chroms(self) -> dict[str, int]:
        return {name: int(length) for name, length in zip(self._af.references, self._af.lengths)}

    def classify(self, read) -> str:
        """'unique' | 'multi' | 'ambiguous' under the configured policy."""
        if self.mapping_method == "nh":
            if read.has_tag("NH"):
                return "unique" if read.get_tag("NH") == 1 else "multi"
            return "ambiguous"
        return "unique" if read.mapping_quality >= self.unique_mapq else "ambiguous"

    def _keep(self, read) -> bool:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            return False
        if self.exclude_duplicates and read.is_duplicate:
            return False
        if read.mapping_quality < self.min_mapq:
            return False
        if self.strand == "+" and read.is_reverse:
            return False
        if self.strand == "-" and not read.is_reverse:
            return False
        if self.mapping_class != "all" and self.classify(read) != self.mapping_class:
            return False
        return True

    def values(self, chrom: str, start: int, end: int) -> np.ndarray:
        out = np.full(end - start, np.nan, dtype="float64")
        chroms = self.chroms()
        if chrom not in chroms:
            return out
        a, b = max(0, start), min(chroms[chrom], end)
        if b <= a:
            return out
        cov = np.zeros(b - a, dtype="float64")
        for read in self._af.fetch(chrom, a, b):
            if not self._keep(read):
                continue
            for blk_start, blk_end in read.get_blocks():
                s = max(blk_start, a)
                e = min(blk_end, b)
                if e > s:
                    cov[s - a : e - a] += 1.0
        out[a - start : b - start] = cov
        return out

    def close(self) -> None:
        self._af.close()


def open_bam_coverage(path: str, name: str | None = None, **kwargs) -> BamCoverageTrack:
    meta_kwargs = {
        k: kwargs.pop(k)
        for k in ("sample", "group", "replicate", "metadata", "display")
        if k in kwargs
    }
    meta = TrackMeta(
        name=name or os.path.basename(path),
        path=os.fspath(path),
        kind=TrackKind.CONTINUOUS,
        format="bam",
        **meta_kwargs,
    )
    return BamCoverageTrack(meta, **kwargs)
