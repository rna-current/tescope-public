# Signal representation: fills, smoothing, distributions, mapping classes

The pipeline is always:

```
RAW SCIENTIFIC VALUES  →  optional visualization transforms  →  graphics
```

Canonical results (`profile_*.tsv`, QC, all analysis tables) are **never**
modified by anything on this page. Display transforms are explicit in the
GUI, the config, and figure provenance (`.figure.json` embeds the full
config, so any smoothing/mapping settings are recoverable per figure).

## Filled profiles

`style.representation: filled` draws each profile line plus a translucent
area to the zero baseline (`fill_alpha`, default fill color = series
color; per-series `fill_color` override available). The plotted values
are identical to `line` mode.

## Smoothing (display-only)

`style.smoothing_method`: `none` (default) | `moving_average` (odd
`smoothing_window`) | `gaussian` (`smoothing_sigma`, in bins) | `savgol`
(odd `smoothing_window` ≥ 5, `smoothing_polyorder`). All NaN-aware: bins
without data stay empty; smoothing never invents values there
(Savitzky–Golay windows containing gaps fall back to the NaN-aware
moving average for that point).

- The **unsmoothed profile remains the canonical numerical result** —
  `profile_*.tsv` is untouched.
- With `smoothing_show_raw: true` (default) the raw profile stays visible
  as a thin low-alpha trace under the smoothed line, so smoothing can
  never quietly hide the observations.
- When smoothing is active, a separately labeled
  `profile_smoothed_*.tsv` is exported carrying the method, parameters,
  source profile and TEscope version.

Do not smooth until peaks appear; the raw trace and the recorded
parameters exist precisely to keep this honest.

## Distribution bands: SEM vs per-copy quantiles

Two different quantities — never confuse them:

- `band: sem` (default): the **uncertainty of the mean** (existing SEM;
  between-replicate SEM for groups).
- `band: iqr`: the **variation across TE copies** — per-bin 25–75%
  quantiles of the per-copy projected values, plus an optional 10–90%
  outer ribbon (`outer_band`). Available for ungrouped continuous tracks
  (a group has no single per-copy matrix; grouped series show no band in
  this mode). Bands use the series color at decreasing alpha and remain
  interpretable in grayscale.
- With `band: iqr`, `profile_quantiles_*.tsv` (p10/p25/p50/p75/p90 +
  n_observations) is exported alongside the canonical profile.

## Copy support

Every figure already carries the aligned **TE copies per bin** panel
(and per-subset support panels in the fragment comparison): apparent
signal supported by few fragments is visibly distinguishable from signal
with deep copy support. Support counts are never mixed into amplitude.

## Mapping classes (BAM/CRAM only)

Per-track **analysis** options (they change the observation set, are
validated up front, and are recorded in config, QC and provenance):

```yaml
tracks:
  - {file: reads.bam, name: unique, mapping_class: unique, mapping_method: nh}
  - {file: reads.bam, name: multi,  mapping_class: multi,  mapping_method: nh}
```

- `mapping_method: nh` — NH tag = 1 → unique; NH > 1 → multi; missing →
  ambiguous. Requires NH tags **in the actual file** (TEscope inspects
  it and fails clearly if absent).
- `mapping_method: mapq` — MAPQ ≥ `unique_mapq` (default 30) → unique;
  below → ambiguous. MAPQ cannot positively identify multi-mappers, so
  `multi` is not available under this method.
- Classification is **aligner-dependent** (NH conventions and MAPQ
  scales differ between aligners); interpret accordingly.
- Secondary/supplementary alignments stay excluded, so a multi-mapping
  read is counted once at its primary position — never multiplied.

When an analysis contains an explicit unique + multi pair from the same
BAM, the figure gains an aligned **fraction unique (of classified
coverage)** panel and a `mapping_support.tsv` export.

**Pre-aggregated tracks (bigWig/bedGraph/BED)**: mapping information does
not exist in the file, and TEscope refuses to pretend otherwise —
mapping options are rejected with a clear message, GUI controls are
disabled with an explanation, and QC records
`mapping policy: unavailable from track`.

## Fragment comparison figure

View mode `fragments` (GUI View selector, `style.fragment_panels: true`
for export) shows **all / full-length / 5′-truncated** copies of the same
analysis on one consensus axis with per-subset copy-support panels. It
regroups the existing per-copy projections by the existing completeness
classes — a subset's profile equals what the corresponding
completeness-filtered analysis produces (tested).

## Heatmap options

`heatmap_sort`: `cons_start` (default) | `genome_length` |
`consensus_coverage` | `total_signal` — row **order** only, values never
change (the title records the sort). `heatmap_show_profile: true` adds
the aggregate profile in an aligned panel above the heatmap.

## Feature-region shading

`annotation_shading: true` washes verified/user feature intervals as a
light background across signal panels (lane colors, low alpha).
Coordinates come verbatim from the feature records — never invented.

## Long-read / Nanopore signals

Modification-probability or long-read coverage tracks are ordinary
continuous inputs (bedGraph/bigWig probabilities, or BAM coverage with
mapping classes where tags allow). A fully synthetic, ground-truth
demonstration ships in `examples/make_nanopore_demo.py`; every file and
figure it produces is labeled **SYNTHETIC NANOPORE DEMONSTRATION — NOT
BIOLOGICAL EVIDENCE**, and that policy applies to any synthetic data:
they demonstrate software behavior, never biology.
