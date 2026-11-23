# Quick start

## GUI

1. **Open TEscope** (`tescope-gui`, or the desktop app once installed).
2. **Choose genome** — e.g. `hg38`. The first time, click *Download
   annotation* (one-time, ~10 min; cached afterwards).
3. **Choose TE** — browse Class → Family → Name, or type `L1HS` in the
   search box. The copy count is shown next to each name.
4. **Drag & drop signal files** — bigWig/bedGraph coverage, BED/narrowPeak
   peak calls, site lists. Assign *Group* and *Rep* in the track table to
   combine replicates (e.g. group `control` reps 1–3, group `treated`
   reps 1–3).
5. **Set filters** (optional) — e.g. *Minimum consensus coverage ≥ 0.5* to
   drop small fragments. The panel shows how many copies remain.
6. **Run.**
7. **Inspect** — Profile tab (stacked or overlay view), Heatmap tab
   (per-copy signal), QC tab (counts, denominators, warnings).
8. **Export** — *Export figure…* (PNG/SVG/PDF), *Export data…* (all result
   tables as TSV), *File → Save project…* (YAML that reproduces the run).

## CLI

```bash
# one-time reference setup
tescope reference fetch hg38

# what's annotated?
tescope te list hg38 --search L1
tescope te info hg38 L1HS

# quick profile
tescope profile \
    --genome hg38 --te L1HS \
    --signal rnaseq_rep1.bw --signal rnaseq_rep2.bw \
    --peaks m6a_peaks.narrowPeak \
    --min-consensus-coverage 0.5 \
    --bin-size 10 \
    --out results/ \
    --save-config l1hs_analysis.yaml

# reproducible re-run (identical to what the GUI executes).
# NOTE: `tescope run <file>` takes a project file that YOU created —
# via --save-config (above), the GUI's File > Save project, or by hand
# (schema below). No analysis.yaml ships with TEscope. For an instant
# self-contained demo, see examples/README.md.
tescope run l1hs_analysis.yaml
```

`results/` then contains per-bin profile tables (`profile_*.tsv`),
occurrence tables, projected peak coordinates, the QC table, the exact
config (`analysis.yaml`), and publication-ready figures with `.figure.json`
metadata sidecars.

## A project file, by hand

```yaml
assembly: hg38
repeat:
  name: L1HS
filters:
  min_consensus_coverage: 0.5
tracks:
  - file: control_rep1.bw
    group: control
    replicate: 1
  - file: control_rep2.bw
    group: control
    replicate: 2
  - file: treated_rep1.bw
    group: treated
    replicate: 1
  - file: m6a_peaks.narrowPeak
    kind: interval
aggregation:
  bin_size: 10
  statistic: mean
output:
  dir: results
```

Track paths are resolved relative to the YAML file's directory.
