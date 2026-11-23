# TEscope User Guide — interpretation concepts

Short, explicit statements of what TEscope does and does not do.
These matter for correct scientific reading of every figure.

## Upstream responsibilities vs TEscope's

TEscope **never** aligns reads, calls peaks, normalizes libraries, or
computes IP/input ratios. You bring processed signal (bigWig/bedGraph
values, BED/narrowPeak intervals, aligned BAM/CRAM); whatever
normalization those files carry (or lack) flows through unchanged.
Comparing two libraries' amplitudes is only meaningful if YOU made
them comparable upstream (e.g. depth normalization). TEscope's job:
project each TE copy's genomic interval onto the repeat-consensus
axis (per-RepeatMasker-hit affine mapping, strand-aware), aggregate
across copies, and show honest denominators.

## Raw values vs display transforms

The canonical per-bin values (per-copy means aggregated across copies)
are computed once per run and never modified by any style, smoothing,
normalization-display, zoom, or theme option — this is pinned by
regression tests. Smoothing and the `Normalize:` combo change the
DRAWN curve only; `Export data…` always includes the raw tables (and
labels smoothed companions separately with their parameters).

## SEM band vs copy-distribution bands

- **SEM** = standard error of the plotted mean (between replicates
  for grouped series): how well the mean is estimated.
- **IQR / 10–90% ribbons** = variation ACROSS TE COPIES: how much
  individual copies differ. Available only for ungrouped continuous
  tracks (needs the per-copy matrix).
These are different quantities; neither is a parametric confidence
interval, and TEscope deliberately does not offer one (typical inputs
have n≤2 replicates, where a CI would imply unsupported precision).
TE copies are genomic repeats, **not** independent biological
replicates — never read a copy ribbon as replicate confidence.

## Missing signal vs zero, and copy support

For bigWig/bedGraph, positions without data are MISSING (excluded
from means — never treated as zero). For BAM coverage, zero reads IS
a measured zero. Each profile row carries `n_elements` (copies whose
genomic span covers the bin) and `n_observations` (copies that
contributed a finite value): the denominators change along the
consensus because truncated copies simply do not have 5′ positions.
Bins supported by a handful of copies are noisy — always read
profiles together with the support panel (the fragments view shows
both per subset).

## Mapping classes (NH/MAPQ) and their limits

Young-TE reads multimap by nature; the aligner upstream already chose
primary placements. TEscope's unique/multi views make the observation
set explicit — `unique (NH=1)`, `multi (NH>1)`, or `unique (MAPQ≥t)`
— but cannot recover true read origins, and the semantics depend on
the aligner that wrote the tags. `multi` requires NH (MAPQ cannot
positively identify multi-mappers). Multi-mapping reads are counted
once, at their primary position; secondary/supplementary records are
never counted.

## Compatibility vs identity (sequence & features)

Dfam or user sequences/features are shown only when the source's
name/alias matches AND its length equals the analysis' consensus
model exactly. That gate prevents coordinate nonsense, but a
name+length match is NOT independent proof of base-level identity —
provenance (accession, release, checksum) is recorded so you can
verify. Incompatible sources are withheld with the reason shown, not
stretched to fit. User-supplied libraries are read locally in place;
their licenses remain the user's responsibility, and TEscope never
bundles or redistributes repeat libraries.

## Coordinates

Internally 0-based half-open; every user-facing readout (cursor,
Go-to, search hits, FASTA headers) is 1-based inclusive. The one
0-based input is the user-feature table (labeled in the panel). The
consensus axis is the RepeatMasker model implied by your annotation —
coordinates on a different consensus model (other releases, other
databases) can be offset; compare positions across studies as
distances with that caveat, not as identities.

## Offline behavior and SVG font portability

After the one-time annotation download (and any Dfam lookups you
request), TEscope runs fully offline from its caches. Exported PDFs
embed the resolved font (genuine Arial where installed, else
Liberation Sans). Exported SVGs keep text as editable text declaring
the Arial-first family chain, so another machine renders them with
its own best match — pixel-identical cross-machine appearance is not
guaranteed; use PDF when that matters.

## Loci, subsets, and upstream estimates (0.2.0)

**A "locus" row is one RepeatMasker annotation record.** TEscope
never merges records; the `rm_id` fragment linkage is displayed so
you can see which records RepeatMasker attributes to one insertion.
A row is therefore not asserted to be a complete biological
insertion.

**Selection is not analysis.** Browsing, sorting, highlighting, and
selecting loci change nothing. Only the explicit "Analyze selected
loci"/"Compare with eligible" actions compute a subset aggregate —
recomputed from the original per-copy values with the same
aggregation as the main profile. The eligible population's aggregate
is never modified, and both denominators stay visible. Subsets you
select interactively are exploratory; TEscope deliberately runs no
automatic significance test on them, because an interactively chosen
subset is not an independent confirmatory finding.

**Upstream per-locus estimates are table context, not signal.** A
count or model estimate per locus (e.g. Telescope's `final_count`)
has no positional resolution: TEscope will sort and compare loci by
it, but never spreads it along the consensus or draws it as a curve.
Declared semantics are enforced: a total that includes unique
evidence is never added to its unique component; a missing component
is *unavailable*, not zero; native scores keep their tool-specific
meaning (no cross-tool confidence score); and a selected subset
never triggers re-estimation or reassignment. Estimates from an
RNA-seq quantifier describe expression evidence under that tool's
model and assumptions — they are not chromatin signal and are
labeled with their assay.

**Feature projection is bounded by the mapping that exists.** The
RepeatMasker `.out` format retains hit endpoints, not per-base
alignments, so consensus→locus projection is affine within each hit:
good enough to show *where* a feature falls on a copy, not to claim
base accuracy. Features outside a copy's footprint get "no verified
correspondence" rather than a stretched interval, and partial
overlaps are flagged. Projection never proves an intact ORF or an
active element; absence from a user-supplied annotation file is
absence from that file, not biological absence.
