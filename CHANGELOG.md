# Changelog

TEscope follows [semantic versioning](https://semver.org).

## [0.2.1] — 2026-11-24

- Fixed: a project opened in the GUI (File > Open project) now keeps the
  file's settings that have no GUI control (aggregation bin size and
  statistic, output settings, chromosome-alias policy, export dpi). In
  0.2.0 such a run used the default bin size (10 bp), so it did not
  reproduce the command-line run of the same file.
- Packaging: the Linux bundle is built on a pinned image (glibc 2.39) and
  the build checks that requirement. The 0.2.0 Linux bundle required
  glibc 2.43 and should not be used.

## [0.2.0] — 2026-11-23

First public release. Adds the Loci tab and Locus Inspector, explicit
locus-subset analysis, user TE feature files (template TSV, BED, GFF3),
upstream per-locus estimate import, schema-2 project files, and
provenance sidecars on toolbar saves. See `RELEASE_NOTES.md`.

## [0.1.0] — development (not released)

- Canonical RepeatMasker parser (native `.out` and UCSC `rmsk` tables)
  preserving repeat-consensus coordinates
- Assembly registry (hg38, hg19, mm39, mm10) with download/cache and
  custom RepeatMasker import
- Fragment-aware consensus-coordinate projection (strand-aware, affine
  within each hit)
- TE filtering: genomic length, consensus coverage, positional completeness
- Continuous tracks (bigWig via pyBigWig/pybigtools, bedGraph), interval
  tracks (BED, narrowPeak, broadPeak), point tracks, optional BAM coverage
- Aggregation with per-bin denominators (n_elements, n_observations),
  group/replicate combining, occurrence profiles for peaks
- QC reporting, chromosome-name reconciliation (chr1/1, chrM/MT)
- Headless CLI (`tescope run/profile/reference/te`) and PySide6 GUI
  sharing one analysis core
- Reproducible YAML project files; TSV/CSV export; PNG/SVG/PDF figures
  with metadata sidecars
- Synthetic golden test dataset, 159 tests, cross-platform CI
- Validated end-to-end against real hg38/hg19/mm39/mm10 UCSC annotations
  and real Nanopore dRNA data (BAM coverage, 10.6M-site bedMethyl points)
- Signal representations: filled profiles, display-only smoothing
  (moving average/Gaussian/Savitzky-Golay, raw+smoothed overlay,
  separately-exported smoothed tables), per-copy quantile ribbons
  (IQR + 10-90%), fragment-comparison figure, heatmap sorting +
  profile panel, feature-region shading, explicit BAM mapping classes
  (NH/MAPQ policies) with fraction-unique panel and honest
  "unavailable from track" for pre-aggregated inputs; input-aware
  capability gating in the GUI; synthetic Nanopore demonstration
- Publication plot-style subsystem: screen/publication presets, curated
  palettes (default/colorblind/grayscale with print-safe line-style
  cycling), per-series colors/line styles/widths, typography, legend,
  grid, display-only axis windows, physical figure sizes, annotation-lane
  and heatmap presentation, SVG/PDF/PNG/TIFF export — a pure rendering
  layer, proven not to change any scientific output
- Interactive consensus browser: wheel zoom / toolbar pan+box-zoom /
  go-to position-range / full view / view history on one synchronized
  consensus axis across all panels; hover readout (1-based base, signal,
  copy support, features) and click position marker — all display-only
- Consensus sequence support: Dfam automatic source (official API,
  CC0, cache-first offline) and user-provided FASTA/library files (read
  locally, never bundled or redistributed; licensed libraries stay under
  the user's own license), gated by the same exact-model compatibility
  verification as feature annotation (verified/unverified/incompatible/
  unavailable; incompatible models such as current split-Dfam L1HS are
  refused, never joined or remapped); zoom-dependent sequence track
  (color strip → vector base letters), copy/FASTA export with
  provenance-rich headers, literal search with zoom-to-hit,
  `tescope sequence show/export/search` CLI, full provenance (accession,
  release, checksum, path) in QC and figure sidecars; sequence controls
  capability-gated with explanatory tooltips
- Polished desktop theme (light default, dark + system via View → Theme):
  centralized token/QSS theme system, original SVG icons, card layout —
  application chrome only; exported figures keep the publication style
- Linux: GUI validated on a real X11 display; headless CLI byte-reproducible;
  PyInstaller bundle with third-party license notices; Windows Inno Setup
  installer script (unsigned; CI-built)
