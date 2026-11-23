# GUI guide

Launch with `tescope-gui` (or `tescope gui`), or the desktop app once
installers ship.

## Layout

- **Header** — application identity and version.
- **Left column** (scrolls on short windows) — the analysis definition:
  reference genome, TE selection, filters, track manager; the **Run
  analysis** button stays pinned below it.
- **Right area** — results: Profile, Heatmap and QC tabs.
- **Status bar** — progress and color-coded state (running / done /
  notes / failed).

## Appearance

*View → Theme* switches between **Light** (default), **Dark**, and
**System** (follows the OS setting where Qt can detect it); the choice is
remembered. The theme styles only the application chrome — exported
scientific figures always use the publication style and are identical
under every theme.

## Reference genome

Pick an assembly from the dropdown. Uninstalled assemblies are labeled
*(not downloaded)*; click **Download annotation** (progress in the status
bar; the download runs in the background and is cached permanently).
**Import custom…** loads your own RepeatMasker `.out` file under a name
you choose — for any organism or assembly.

## Choosing a TE

The tree shows Class → Family → Name with copy counts; the search box
filters by name (e.g. `L1HS`, `AluY`, `MER41`). Select a leaf (a repeat
name / subfamily).

## Filters

Enable the checkboxes you need; the counts below update live:

```text
All copies                 1,523
Length >= 1 kb               812
Consensus coverage >= 50%    391
```

None enabled = analyze every copy. See [filtering.md](filtering.md) for
what each filter means.

## Tracks

Drag & drop files anywhere on the window, or click **Add…**. Per row you
can edit:

- **Name** — display name (double-click)
- **Group** — tracks with the same group are combined as replicates
- **Rep** — replicate number
- **Show** — include/exclude from plots without removing

Formats are detected from the extension; unsupported files produce an
explicit message, not silence. Reorder with Up/Down; remove with Remove.

## Consensus features

The *Consensus features* section annotates the profile's consensus axis.
Enable it, and TEscope looks up the selected repeat in Dfam (cache first;
*Check Dfam…* fetches on demand) — the status line states whether a
record is verified, incompatible with your consensus, or absent. Add your
own features in the table (0-based half-open coordinates) or *Import…*
a TSV/YAML file; they always appear in a separate "user-defined" lane.
See docs/annotations.md — TEscope never invents annotations.

## Consensus sequence

The *Consensus sequence* section chooses the sequence source: Automatic
(Dfam, cache-first), a user-provided library or FASTA (read locally in
place — never uploaded or copied), or none. *Check source…* verifies
compatibility with the analysis' consensus model before any base is
shown; incompatible or unavailable sources are explained in the status
line and leave the analysis untouched. See docs/sequence-browser.md.

## Running and results

**Run analysis** executes in a background thread (the window stays
responsive; progress in the status bar). Errors appear as dialogs with
actionable messages ("Reference annotation not yet downloaded", "Track X
shares no chromosomes with the selected assembly", "Only 7 TE copies
remain after filtering" as a note).

- **Profile tab** — stacked (one panel per group/track) or overlay view;
  y-axis independent/shared; optional per-profile normalization (z-score /
  min–max — recorded in the project file); SEM band; individual
  replicates. The bottom panel always shows TE copies per bin. The view
  is an interactive consensus browser: wheel zoom, toolbar pan/box-zoom/
  history, *Go to* a base or range, hover readout, click marker, and the
  sequence controls (show/copy/export/search) when a verified consensus
  sequence exists — see docs/sequence-browser.md.
- **Heatmap tab** — copies × bins matrix per continuous track.
- **QC tab** — the full QC table (also exportable as TSV via *Export
  data…*).

## Export & projects

- **Style…** — the figure-style dialog: preset (screen/publication),
  palette, per-series colors (picker or hex) and line styles/widths,
  typography, legend, grid, display ranges, annotation-lane and heatmap
  appearance. Apply redraws live; everything saves into the project file
  and reproduces headlessly. Rendering only — results never change.
- **Export figure…** — PNG, SVG, PDF or TIFF, with a `.figure.json`
  reproducibility sidecar.
- **Export data…** — writes every result table plus `qc.tsv` and
  `analysis.yaml` to a chosen directory.
- **File → Save project… / Open project…** — the YAML project file
  restores genome, TE, filters, tracks and plot options, and runs
  identically through the CLI: `tescope run analysis.yaml`.
