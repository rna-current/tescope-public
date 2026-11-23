# Output files and interpretation

Every analysis produces machine-readable tables sufficient to reconstruct
the figures independently — no result exists only as an image.

## Tables

### `profile_<track>.tsv` — per-bin continuous profile

| column | meaning |
|---|---|
| `repeat_name`, `signal_name`, `group`, `aggregation` | context |
| `consensus_bin` | bin index |
| `consensus_start`, `consensus_end` | bin interval on the consensus (bp, 0-based half-open) |
| `mean`, `median` | across copies of the per-copy bin means |
| `sem` | standard error across copies (NaN when < 2 observations) |
| `n_elements` | copies whose consensus footprint covers the bin |
| `n_observations` | copies that contributed ≥ 1 finite base value |

**Always read `n_elements`.** In fragmented families the left half of the
profile may aggregate 90 copies and the right half 900; identical y-values
do not mean identical evidence.

### `group_<group>.tsv`

Replicate-combined profile: mean of the replicates' statistic per bin,
between-replicate `sem`, `n_replicates`, and conservative (minimum)
denominators.

### `occurrence_<track>.tsv` — peaks per bin

- `peak_count`: projected (copy, peak) intervals overlapping the bin. One
  genomic peak overlapping 3 copies counts 3 times — the unit is the
  copy-projection.
- `copies_with_peak` and `fraction_copies_with_peak`: how many / what
  fraction of the bin-covering copies have ≥ 1 peak there. This is the
  "what fraction of L1HS copies carry this mark at this position" number.

### `projected_intervals_<track>.tsv`

One row per (TE copy, overlapping peak): original genomic peak
coordinates, projected consensus coordinates (clipped to the copy), name,
score, strand — the raw material for custom plots.

### `points_<track>.tsv`, `point_occurrence_<track>.tsv`

Projected single-base events with fractional consensus positions, plus
per-bin counts and per-element rates.

### `copies.tsv`

The post-filter copy set with all canonical and derived columns
(consensus coverage, completeness class, ...). Use it to audit exactly
which elements your profile aggregates.

### `qc.tsv`

The full QC summary: filter step counts, per-track chromosome
reconciliation, signal missingness, peak counts, warnings.

### `analysis.yaml`

The exact configuration that produced this directory, including annotation
provenance (source URL/path, retrieval date, MD5) and TEscope version.
`tescope run analysis.yaml` reproduces the analysis.

## Figures

- `profile_stacked.*` — one panel per group/track, plus interval and point
  panels, plus the **TE copies per bin** panel at the bottom (the
  denominator, always visible).
- `profile_overlay.*` — continuous signals overlaid on one axis (mind the
  normalization setting).
- `heatmap_<track>.*` — copies × bins matrix, copies sorted by consensus
  start; the most honest view of fragment structure.

Each figure is written in the configured formats (PNG/SVG/PDF) with a
`.figure.json` sidecar recording the version, config and provenance needed
to regenerate it.
