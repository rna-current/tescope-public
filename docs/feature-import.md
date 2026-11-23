# Importing TE features and upstream locus estimates

Two optional input kinds extend an analysis. Both are display/
context layers or locus-table columns — neither changes the computed
consensus profile.

## 1. TE features / annotations (TSV template)

One file = one coordinate target, declared in a header comment:

```
# tescope-features v1
# target: consensus            (or: genome)
# reference: L1HS              (consensus: model name; genome: assembly e.g. hg38)
# coordinates: 0-based-half-open   (TSV template is always this)
feature_id	start	end	strand	label	type	parent_id	source
ORF1	1987	3004	+	ORF1	CDS		L1 literature
ORF2	3067	6889	+	ORF2	CDS		L1 literature
5UTR	0	909	+	5'UTR	UTR		L1 literature
```

- `feature_id` unique per file; `parent_id` optional (must reference
  an existing feature_id).
- `target: consensus` features attach to the consensus axis (same as
  the GUI user-feature table). `target: genome` features attach to
  genomic coordinates and appear in the Locus Inspector for loci they
  overlap (assembly must match the analysis assembly — no silent
  aliasing).
- Coordinates are nucleotides. Amino-acid coordinates are not
  supported and are rejected with an explanation (convert them
  yourself only if you know the CDS frame; TEscope does not guess).
- BED (0-based half-open) and GFF3 (1-based inclusive) are also
  accepted for `target: genome`; the convention is taken from the
  format specification, never inferred from values.

Validation rejects: unknown/duplicate IDs, end <= start, negative
coordinates, out-of-bounds vs a known reference length, unknown
strand values, dangling parent_id, contig names absent from the
assembly, amino-acid units.

## 2. Upstream per-locus estimates (TSV)

For locus-level TABLE context only (sorting, comparison). Estimates
are never drawn as positional signal: a per-locus count cannot
produce a coverage profile.

```
# tescope-locus-estimates v1
# tool: telescope 1.0.3        (free text, recorded verbatim)
# assay: RNA-seq               (required; shown in the UI)
# units: fragments             (required; shown in the UI)
# total_includes_unique: true  (required: true/false/unknown)
# assembly: hg38
locus_id	total	unique	score
chr1:1000-7000(+)	120	45	0.97
```

- `locus_id` must be a TEscope stable locus ID
  (`chrom:start-end(strand)`, 0-based half-open as exported by the
  Loci panel / copies table) in the declared assembly; unmatched IDs
  are reported, not dropped silently.
- `total` and `unique` keep native units. When
  `total_includes_unique: true` (e.g. Telescope `final_count` with
  `unique_count`), TEscope never adds them. A missing `unique`
  column is *unavailable*, not zero. `score` keeps its native
  meaning (recorded from the header), no cross-tool normalization.
- Native Telescope reports (`telescope-stats_report.tsv`) can be
  converted with the documented column mapping (transcript→locus_id
  via your annotation; final_count→total; unique_count→unique;
  final_conf→score). TEscope does not run or re-fit the upstream
  model, and selecting a locus subset never triggers reassignment.
