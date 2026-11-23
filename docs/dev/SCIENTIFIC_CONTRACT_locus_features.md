# Scientific contract — locus exploration, user TE features, upstream estimates

Governs the 2026-10-08 extension (dev/locus-features). UI behavior is
wired to these rules; anything not derivable from them is omitted
from the UI rather than approximated.

## 1. Entities and identity

- **Annotation record (fragment)**: one RepeatMasker hit row.
  Identity = (assembly, chrom, start, end, strand, repeat name,
  annotation source + retrieval metadata). Records are never
  silently merged or deduplicated; grouping of fragments into a
  larger "copy" exists only where the already-validated per-copy
  grouping logic produces it, and the fragment IDs remain visible.
- **Locus** (as shown in the Loci panel) = one analyzed copy after
  the existing grouping — NOT asserted to be a complete biological
  insertion. The panel shows fragment count and completeness class so
  the user can see what the row actually is.
- **Stable locus ID**: `chrom:start-end(strand)` in the annotation's
  assembly, 1-based inclusive for display, stored internally 0-based
  half-open — stable across sessions for the same annotation input.
  Selection export/import uses this ID + assembly + annotation
  provenance; a selection fails loudly (listed, not dropped silently)
  if an ID no longer matches the current eligible population.

## 2. Populations and state

- **Eligible population**: copies passing the scientific filters of
  the current run. This is the denominator of every aggregate unless
  an explicit subset analysis is applied.
- **Highlighted locus**: browsing focus; zero analytic effect.
- **Selected loci**: a candidate subset; zero analytic effect until
  the explicit action.
- **Applied subset**: result of "Analyze selected loci" /
  "Compare subset vs eligible" — a separate, labeled computation
  over original per-locus values; the eligible population remains
  visible as context. No upstream estimate is re-fit or reassigned
  against a subset, ever.
- Sorting, clicking, zooming, styling, smoothing: display-only.

## 3. Signal classes (never conflated)

| class | example | positional? | locus-of-origin confidence |
|---|---|---|---|
| Supplied continuous signal ("As supplied") | ENCODE fold-change bigWig | yes | whatever the upstream pipeline did — NOT re-classified by TEscope |
| Alignment-derived coverage with mapping policy | BAM, NH=1 / NH>1 / MAPQ>=t | yes | conditional on aligner, reference, tags; policy shown, never a biological guarantee |
| Interval/peak evidence | narrowPeak | interval | as supplied |
| Upstream locus estimate | Telescope `final_count` per locus | **no** | per the upstream method's own semantics |

Rules:
- BigWig/bedGraph values are never labeled unique/multi by TEscope;
  only user-supplied provenance text can describe them.
- Upstream per-locus estimates support tables, sorting, and
  locus-level comparison ONLY. They are never spread over positions,
  never drawn as a profile, and never combined on an axis with
  coverage units.
- Telescope semantics (verified against official docs 2026-10-08):
  `final_count` is a total that already includes unique evidence —
  `unique_count` is a component, never an addend. Missing
  decomposition = unavailable, not zero. Native scores
  (`final_conf`, `final_prop`) keep their method-specific meaning;
  no common confidence score is manufactured.
- Assay binding: an estimate table declares its assay (e.g. RNA-seq);
  the UI labels it and does not present it as chromatin signal.

## 4. Coordinates

- Internal: 0-based half-open everywhere. Display: 1-based inclusive.
- Two feature targets: (A) consensus/reference-model features;
  (B) genomic (assembly) features. Each imported file declares its
  target, reference identity, and coordinate convention (BED =
  0-based half-open; GFF3 = 1-based inclusive; TSV template declares
  explicitly). Conventions come from the format, never guessed from
  values.
- Consensus→locus projection of a feature is offered only through the
  existing verified per-fragment affine mapping, only over the
  fragment intervals the mapping actually covers, and is labeled
  "projected"; where a feature falls outside mapped intervals, the
  locus view shows "no verified correspondence" instead of a drawn
  interval. RepeatMasker endpoints do not justify gap-aware
  base-accurate claims, and the projection is labeled with that
  limitation.
- Amino-acid coordinates are rejected with an explanation (no ×3
  conversion).

## 5. Unknowns and compatibility

- Unknown metadata (aligner, normalization, units, assay) is stored
  and displayed as unknown — never defaulted to a plausible value.
- name+length compatibility (existing sequence gate) is labeled
  compatibility, not identity. Checksums are compared only when both
  sides supply one; an isolated checksum proves nothing and is shown
  as "recorded, not compared".
- Assembly/contig mismatches in feature or estimate imports fail
  validation; aliasing requires an explicit user mapping.

## 6. Aggregation and comparison

- Per-bin aggregate = mean of per-copy means (existing, pinned by
  golden tests); support = n_elements/n_observations, unchanged.
- Feature/domain summaries exist only where positional signal AND
  justified correspondence exist; they document weighting
  (per-locus, equal weight over observed bases), units (input
  units), effective observed length, and included/partial/excluded
  locus counts. Structural absence ≠ 0; missing data ≠ 0.
- Overlapping/nested features are summarized per feature; no derived
  "non-overlapping total" is produced.
- Subset vs eligible comparison recomputes from original per-locus
  values with identical binning/units; denominators of both sides are
  displayed. No automatic significance testing on interactive
  subsets; SEM/IQR ribbons are never called confidence intervals.
- Copy-to-copy variation is labeled as such and kept distinct from
  replicate variation.

## 7. Provenance

Everything needed to recreate an analysis state persists in the
project file (schema-versioned): sources + checksums where computed,
reference/annotation identity, mapping policy, estimate-table
semantics as imported, filters, selected locus IDs, feature inputs,
applied subset, aggregation, and figure-relevant view config. Loading
an rc2 project remains valid; absent new fields mean "feature not
used", never a reinterpretation.
