# TE copy filtering

Genomic TE annotations are dominated by fragments. TEscope's filters let
you choose the copies your question actually concerns, and always show the
counts before/after each step:

```text
L1HS annotations
All copies                    1,523
Length >= 1 kb                  812
Consensus coverage >= 50%       391
Completeness: full length        87
```

## Filters

### Minimum genomic length

`min_genome_length` (bp): drops hits whose genomic footprint is shorter
than the threshold.

### Minimum consensus coverage

`min_consensus_coverage` (0–1): the fraction of the repeat model the hit
spans —

```text
(cons_end - cons_start) / consensus model length
```

A 300 bp hit of a 6 kb LINE-1 model has coverage 0.05 regardless of its
genomic length.

### Positional completeness

Classifies each hit by which consensus ends it reaches, within a tolerance
(default: 5 % of the model length or 20 bp, whichever is larger —
configurable via `end_tolerance_frac` / `end_tolerance_bp`):

| Class | Reaches 5' end | Reaches 3' end |
|---|---|---|
| `full_length` | yes | yes |
| `five_prime_truncated` | no | yes |
| `three_prime_truncated` | yes | no |
| `internal_fragment` | no | no |

Select any subset (e.g. keep only `full_length` for promoter-focused
questions, or only `five_prime_truncated` to study L1 3'-end fragments).
The classification is purely geometric and family-agnostic — no biology is
hard-coded; tune tolerances per analysis when a family needs different
behavior.

### Other options

- `strands: ["+"]` / `["-"]` — restrict genomic strand (rarely needed;
  projection is already strand-aware).
- `include_secondary: false` — drop hits RepeatMasker marked `*`
  (overlapped by a higher-scoring hit). They are kept by default and
  flagged in the `is_secondary` column.

## Interaction with denominators

Filtering changes which copies exist; the per-bin denominators
(`n_elements`, `n_observations`) always describe the **post-filter** copy
set. After aggressive filtering few copies may remain — TEscope warns when
fewer than 10 copies survive ("Only 7 TE copies remain after filtering"),
because aggregate profiles over tiny copy sets are noisy.
