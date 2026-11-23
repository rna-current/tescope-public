# TEscope examples

## Demo dataset (no real data required)

`make_example_data.py` generates a deterministic synthetic dataset that
mimics an L1-like family: a 6 kb consensus, 60 copies on both strands with
realistic fragment structure (mostly 5'-truncated), two RNA-seq-like
conditions (2 replicates each; the "treatment" adds a 5' signal bump), a
broad H3K9me3-like track enriched on full-length copies, narrowPeak calls
over consensus positions 800–1000, and point "modification sites"
clustered at 5.5–6 kb.

```bash
python examples/make_example_data.py example-data

tescope reference add --assembly demo \
    --repeatmasker example-data/demo.fa.out \
    --chrom-sizes example-data/demo.chrom.sizes

tescope run example-data/analysis.yaml
# -> demo-results/ : profile tables, occurrence tables, QC, figures
```

Or explore the same files in the GUI: `tescope-gui`, choose assembly
`demo`, select `L1demo`, drag the bedGraph/narrowPeak/bed files in, Run.

What to look for in `demo-results/profile_stacked.png`:

- **control** group: smooth 3'-biased coverage ramp
- **treated** group: the same ramp plus the designed bump at ~900 bp
- **H3K9me3**: higher over the 5' half — bins there are covered mainly by
  full-length copies, which carry the enrichment
- **5' peaks**: occurrence step at consensus 800–1000 with the projected
  peak rug above it
- **modification sites**: histogram clustered near 5.5–6 kb
- **TE copies per bin**: rises toward the 3' end — the fragment structure
  every fragment-unaware plot would hide

## Real-data recipe

Any coverage bigWig + peak calls produced on hg38/hg19/mm39/mm10 work
directly:

```bash
tescope reference fetch hg38
tescope profile --genome hg38 --te L1HS \
    --signal your_coverage.bw --peaks your_peaks.narrowPeak \
    --min-consensus-coverage 0.5 --out results/
```

Read `docs/limitations.md` (especially multi-mapping) before interpreting
profiles of young TE families.
