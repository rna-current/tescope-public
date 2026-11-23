# Scientific limitations

Read this before publishing TEscope-derived figures. The GUI's simplicity
does not remove any of these caveats.

## 1. Multi-mapping reads (the big one)

Short reads from young, high-identity TE families (L1HS, AluY, recent
ERVs) frequently align equally well to many copies. Standard pipelines
either discard such reads or assign them arbitrarily/randomly. Any
coverage-type signal you load — bigWig from RNA-seq, ChIP-seq, ATAC-seq —
already carries that upstream decision, and TEscope **does not correct
it**: profiles over young families under-represent (MAPQ-filtered) or
smear (random-assignment) the true per-copy signal. TEscope is a
visualization/aggregation tool, **not** a TE expression quantifier; for
locus-level expression claims use dedicated methods (e.g. EM-based TE
quantifiers) and long reads where possible. BAM-derived coverage in
TEscope obeys your MAPQ filter and states this caveat in QC output.

## 2. Annotation dependence

Everything is conditional on the RepeatMasker annotation you load: masking
sensitivity, library version (Repbase vs Dfam eras), and assembly build
all change copy sets and consensus coordinates. TEscope records annotation
provenance (source, date, MD5) in every output; compare analyses only
within one annotation.

## 3. Fragmented and composite annotations

RepeatMasker often splits one biological insertion into several hits
(nested insertions, large deletions) and joins them only via the ID
column. TEscope treats each hit independently; a copy interrupted by an
Alu appears as two fragments. Completeness classes describe hits, not
reconstructed insertions.

## 4. Within-copy linear interpolation

`.out` files lack per-base alignments, so consensus positions within a hit
are interpolated linearly (indels distributed proportionally). At default
bin sizes (≥ 10 bp) the error is usually below one bin for typical
divergences; do not over-interpret single-bin features in highly diverged
families (consider larger bins there).

## 5. Consensus model length is annotation-derived

Model length = max(`end + left`) across the repeat's hits. If a
RepeatMasker run mixed library versions for one name, coordinates may be
slightly inconsistent between hits; TEscope uses the maximum and clips
coverage at 1.0.

## 6. Aggregate profiles are population summaries

A mean profile over hundreds of copies describes the family average in
bulk data — not any specific locus, cell type or allele. Per-bin
denominators (`n_elements`) differ along the consensus in fragmented
families; regions covered by few copies are noisier. The bottom
"TE copies per bin" panel and the heatmap exist so this is never hidden.

## 7. Signal comparability

TEscope never normalizes across tracks silently. If your bigWigs are not
comparably normalized (CPM/RPGC/spike-in), between-track and between-group
comparisons inherit that. The display normalizations (z-score, min–max)
are per-profile visual aids, recorded in the config, and change the y-axis
meaning.

## 8. Assembly mismatch

Signals produced on a different assembly than the selected annotation are
wrong even when chromosome names match. TEscope can only detect naming
mismatches, not liftover errors — it refuses tracks with alien chromosome
names but cannot see a hg19 bigWig labeled with hg38-style names.

## 9. GUI validation status

The full GUI workflow (reference selection, TE search, live filter
counts, track loading, threaded analysis on real hg38 data, results/
heatmap/QC tabs, export, project save/reload) has been validated on a
real Linux X11 display. Windows and macOS GUI behavior is exercised only
through CI's offscreen Qt platform, which validates wiring, not real
usability; real-device validation on physical Windows and macOS machines
is still required before claiming polished desktop support there.

## 10. Consensus feature annotations are display metadata

Database features (Dfam) are shown only for consensus models whose length
exactly matches the source record; homology evidence keeps its source
terminology and is not a gene/ORF claim; many families have no database
features at all, and recent Dfam releases split some classic families so
their features are honestly unavailable for monolithic RepeatMasker
models. User-defined features are the user's responsibility. Annotations
never alter numerical results (docs/annotations.md).

## 11. Mapping classification is aligner-dependent

Unique/multi classification (BAM/CRAM only) relies on NH tags or MAPQ
conventions that differ between aligners; TEscope records the exact
policy but cannot standardize aligner semantics. For pre-aggregated
tracks the mapping policy is genuinely unknowable and reported as such.

## 12. Point tracks in the GUI

Point tracks are fully supported in the core and CLI (`--points`); the GUI
track manager currently loads BED-like files as interval tracks by
default and does not yet offer a "load as points" toggle.

## 13. Consensus sequence display is provenance-gated

TEscope displays a consensus sequence only when its compatibility with
the analysis' RepeatMasker-derived coordinate model has been positively
verified (exact length). This is deliberately conservative: two models of
the same length are not guaranteed base-identical if they come from
different curation lineages — the provenance (accession, release,
checksum, or your own library file) is part of the claim and is recorded
everywhere the sequence is used. For families whose current Dfam model
differs from the annotation model (e.g. L1HS with UCSC hg38 annotation),
the honest state is "incompatible": supply the exact consensus library
your annotation was produced with to enable sequence display. Length
equality also cannot detect a same-length but different-sequence
substitution; the recorded checksum and source identity exist so you can
audit exactly which sequence was shown.
