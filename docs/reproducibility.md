# Reproducibility

## The project file is the analysis

Every TEscope analysis — GUI or CLI — is fully described by one YAML file.
The GUI's *Save project* and the CLI's `--save-config` write it; `tescope
run` executes it; *Export data* includes the exact `analysis.yaml` that
produced the results. GUI and CLI share one implementation
(`tescope.core.analysis.run_analysis`), so a saved GUI session reruns
bit-identically headless.

## Schema

```yaml
tescope_version: 0.1.0          # stamped on save
assembly: hg38
repeat:
  name: L1HS
filters:                        # all optional
  min_genome_length: 1000
  min_consensus_coverage: 0.5
  completeness: [full_length]   # subset of the four classes
  end_tolerance_frac: 0.05
  end_tolerance_bp: 20
  include_secondary: true
  strands: []                   # [] = both
tracks:
  - file: control_rep1.bw       # relative to this YAML's directory
    name: RNA-seq control 1
    kind: continuous            # continuous | interval | point ('' = sniff)
    group: control
    replicate: 1
    visible: true
    color: ""                   # display only
  - file: peaks.narrowPeak
    kind: interval
aggregation:
  bin_size: 10                  # consensus bp per bin
  statistic: mean               # mean | median
plot:
  show_sem: true
  y_axis: independent           # independent | shared
  normalization: none           # none | zscore | minmax (display only)
  show_replicates: false
  dpi: 300
output:
  dir: results
  table_formats: [tsv]          # tsv | csv
  figure_formats: [png, svg]    # png | svg | pdf
chrom_alias: true               # allow chr1<->1 / chrM<->MT reconciliation
annotation_provenance: {...}    # filled in at run time (source, date, md5)
```

Unknown keys are rejected with the offending key named — a typo cannot
silently change an analysis.

## What is recorded where

| Artifact | Records |
|---|---|
| `analysis.yaml` | full config + TEscope version + annotation provenance |
| `qc.tsv` | filter counts, per-track reconciliation, missingness, warnings |
| `*.figure.json` | version, config, provenance, copy counts per figure |
| `provenance.json` (cache) | annotation source URL/path, date, MD5 |

## Numerical reproducibility

- The regression snapshot (`tests/golden/expected_profile.tsv`) pins the
  golden-dataset numbers; CI verifies them on Linux, Windows and macOS on
  every commit (validation milestone I).
- Both bigWig backends (pyBigWig, pybigtools) are tested to produce
  identical values, so platform-dependent backend choice does not change
  results.
- Aggregations are deterministic (no sampling, no seeds).

## Reproducing a published figure

1. Obtain the `analysis.yaml` (or `.figure.json`) alongside the figure.
2. Install the recorded TEscope version: `pip install tescope==<version>`.
3. Verify the annotation MD5 matches your cached copy
   (`provenance.json`), or `tescope reference add` the recorded source
   file.
4. `tescope run analysis.yaml`.
