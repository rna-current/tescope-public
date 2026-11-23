# Input formats

## Continuous / quantitative tracks

Per-base quantitative signal, projected and averaged into consensus bins.

| Format | Extensions | Notes |
|---|---|---|
| bigWig | `.bw`, `.bigWig` | Preferred for genome-wide signal (indexed random access). Read via pyBigWig, or pybigtools on Windows — results are identical (tested). |
| bedGraph | `.bedGraph`, `.bg`, `.bdg` (+`.gz`) | Loaded fully into memory; fine for typical files, use bigWig for very dense genome-wide data. Overlapping intervals are rejected. |
| BAM/CRAM | `.bam`, `.cram` | **Convenience only** (`tescope[bam]`, Linux/macOS). Simple filtered per-base coverage; see caveats below. |

Typical contents: RNA-seq coverage, ChIP-seq / CUT&Tag / CUT&RUN signal,
ATAC-seq signal, m6A/m5C enrichment ratios, Nanopore-derived quantitative
scores, methylation fractions.

Missing data are preserved as missing (NaN), never treated as zero. Bins a
copy covers without data count toward `n_elements` but not
`n_observations`.

## Interval / peak tracks

Kept as intervals — **never converted into curves**.

| Format | Min columns | Score used |
|---|---|---|
| BED | 3 | column 5 (if numeric) |
| narrowPeak | 10 | signalValue (column 7) |
| broadPeak | 9 | signalValue (column 7) |

Each peak overlapping a TE copy is projected to that copy's consensus
coordinates (clipped to the copy). Outputs distinguish **peak count per
bin** from **fraction of covering copies with a peak** — the denominators
differ and both are reported.

## Point tracks

Single-base events: RNA modification sites, A-to-I editing sites, CLIP
crosslink sites, per-base Nanopore modification calls. Load any BED-like
file with `kind: point` (CLI `--points`, GUI: not yet exposed — see
limitations). 1 bp intervals are used as-is; longer intervals contribute
their midpoint. The score column is carried through as the event value.

## BAM/CRAM caveats

- Input must be coordinate-sorted and indexed; a missing index is reported
  (`BAM index missing`) and can be created automatically when enabled.
- MAPQ filtering, duplicate filtering and strand selection are exposed.
- **Multi-mapping reads are not rescued.** Coverage inside young, high-copy
  TE families under-represents true signal when reads cannot be uniquely
  placed. TEscope does not quantify TE expression — see
  [limitations](limitations.md).

## Chromosome naming

Track chromosome names are reconciled with the reference: exact matches
first, then unambiguous aliases (`chr1`↔`1`, `chrM`↔`MT`↔`M`). Anything
ambiguous is reported, never guessed. A track sharing no chromosomes with
the selected assembly stops the analysis with an explicit error naming the
file. Alias reconciliation can be disabled (`chrom_alias: false`).
