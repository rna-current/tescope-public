# TEscope — offline quick help

TEscope visualizes and aggregates your genomic signal files (bigWig,
bedGraph, BED/narrowPeak, BAM coverage) on transposable-element
consensus coordinates, fragment- and strand-aware. Your data stay on
this machine; network access is used only to download public reference
annotation and — when you request them — public Dfam consensus
features/sequences. Cached references keep working offline.

## Quick start (GUI)
1. Reference genome: pick an assembly (hg38/hg19/mm39/mm10) and press
   "Download annotation" once (network required), or "Import custom…"
   to use your own RepeatMasker .out / UCSC rmsk table.
2. Transposable element: search and select a repeat (e.g. L1HS).
3. Optional TE copy filters: length / consensus coverage /
   completeness classes.
4. Add signal tracks (File ▸ Add track… or the track panel), assign
   groups/replicates if desired.
5. Run analysis. Explore Profile / Heatmap / QC tabs; View =
   overlay / stacked / fragments; Style… for presets and bands.
6. Export figure… (PNG/SVG/PDF/TIFF), Export data… (TSV), and save
   the configuration (File ▸ Save config) for a reproducible rerun
   (CLI: `tescope run config.yaml`).

## Key concepts
- Bin values are per-copy means; `n_elements` / `n_observations`
  report honest per-bin copy support.
- SEM bands = uncertainty of the mean; IQR/decile bands = variation
  across TE copies — not interchangeable.
- Smoothing is display-only; exported tables stay raw.
- BAM mapping classes (unique/multi) follow NH or MAPQ policies and
  are gated by what the file actually supports.

## More
Full documentation and the User Guide ship with the source repository
(docs/) and the project website when available. License: see the
License button (MIT). Citation: see "How to cite".
