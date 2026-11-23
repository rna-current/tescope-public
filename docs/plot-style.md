# Plot style and figure customization

TEscope's plot style is a **pure rendering layer**: every option below
changes only how figures look. No style setting can alter projected
coordinates, binned values, aggregation, `n_elements`/`n_observations`,
normalization, filtering, feature coordinates or any exported scientific
table — this is pinned by regression tests that render the same analysis
under heavily different styles and compare all result tables.

## Presets

| preset | intent |
|---|---|
| `screen` (default) | comfortable GUI viewing; identical to TEscope's established look — old configs render unchanged |
| `publication` | manuscript-oriented: ~17.8 cm (double-column) width, restrained typography, thinner lines |
| `publication_minimal` | like publication with quieter typography and thinner lines, for dense multi-panel layouts |
| `presentation` | slides/projector: 24 cm width, large type, thick lines; *Reset* also selects the filled representation |
| `high_contrast` | screen width with enlarged type and thick (2.4 pt) lines for low-vision/projector readability |
| `custom` | same base as screen; signals user overrides |

Uncertainty display remains `band: sem` (uncertainty of the mean) or `band: iqr` (per-copy variation, with the 10–90% outer ribbon). A parametric confidence-interval band is deliberately NOT offered: typical TEscope inputs carry n≤2 replicates, where a CI would imply precision the data cannot support — the honest alternatives above are shown instead.

A preset sets defaults; any individual option can be overridden. *Reset
to preset* (GUI) or a fresh `style:` section restores defaults. Figure
widths of ~8.5 cm / ~17.8 cm approximate common single/double journal
columns — no specific journal standard is claimed.

## Series colors, line styles, widths

Each plotted series (a group, or an ungrouped track) is keyed by its
**stable name** — never by row number — and can carry:

- `color`: any hex or named matplotlib color (validated; invalid values
  are rejected with the series named). Explicitly chosen colors are
  **never reassigned** between runs.
- `line_style`: `solid` / `dashed` / `dotted` / `dashdot`
- `line_width`, `alpha`

Unstyled series draw from the selected palette in display order:
`default` (Okabe-Ito, the established TEscope cycle), `colorblind`
(Paul Tol bright), or `grayscale` — which automatically cycles line
styles so gray series stay distinguishable in print. (Palettes are
chosen for color-vision friendliness; no accessibility certification is
claimed.) Stale style entries for series that no longer exist are kept
harmlessly and simply unused.

The SEM uncertainty band keeps its series' color; its opacity is
`band_alpha`. Band visibility remains the existing `plot.show_sem`
switch; the statistics behind it are untouched.

## Typography, titles, legend, axes

All figure text uses the declared family chain **Arial → Liberation
Sans → DejaVu Sans → installed CJK fallbacks → sans-serif**,
identically in the GUI preview, CLI exports, and every style preset.
Hosts without native Arial (typical Linux) render with Liberation Sans
(metrically Arial-compatible) — the declared preference stays Arial,
and SVG exports carry the Arial-first `font-family` metadata (text
remains editable text, `svg.fonttype "none"`). PDF embeds the
*resolved* face as TrueType (`pdf.fonttype 42`): real Arial where
installed, otherwise the first available fallback — a backend
property, documented rather than worked around by rasterizing.

- Font sizes (pt; 0 = preset default): `title_size`, `label_size`,
  `tick_size`, `legend_size`, `annotation_size`. System fonts only —
  nothing is bundled; missing families fall back predictably (DejaVu
  Sans primary, host-installed CJK-capable fallbacks for non-Latin
  labels).
- `title`, `x_label`, `y_label` override presentation text; defaults
  keep TEscope's semantics (repeat/assembly/copy counts). TEscope never
  invents biological titles.
- Legend (overlay view): `show_legend`, `legend_position` ∈ best /
  upper right (default) / upper left / lower right / lower left /
  outside right. The outside position exports with a tight bounding box
  so it is never clipped.
- `show_grid` + `grid_alpha`.
- **Display limits** `y_min`/`y_max` (signal panels) and `x_min`/`x_max`
  (consensus view window, bp) are *view windows only*: the underlying
  data and all exported tables keep the full range. Analytical
  restriction is what the `filters:` section does — the two are
  deliberately separate.

## Figure size

`fig_width_cm` / `fig_height_cm` in centimetres; 0 means automatic
(width from preset; height scales with panel count, annotation lanes
included). Saved dimensions reproduce exactly.

## Annotation lane style

`annotation_alpha`, `annotation_border_width`, `annotation_size`, and
lane colors (`database_lane_color`, `user_lane_color`). Origin always
stays visible: the user-defined lane keeps its hatch and both lanes keep
their labels regardless of chosen colors, and styling can never move a
feature — coordinates come from the feature records alone.

## Heatmap

`heatmap_cmap` from a curated list (viridis, magma, cividis, Greys,
RdBu_r, coolwarm), optional `heatmap_vmin`/`heatmap_vmax`, and
`heatmap_colorbar`. **Color limits are display-only — they are not
normalization**; the underlying per-copy matrix is never rescaled.

## Export

Formats: **SVG and PDF (true vector — never silently rasterized)**,
PNG and TIFF (raster; DPI applies, e.g. 300/600/custom via `plot.dpi`).
DPI is meaningless for vector output and ignored there. Convenience
combinations: vector = `[pdf, svg]`; high-res raster = `[png, tiff]` at
600 dpi; screen = `[png]` at 150 dpi.

## Signal representation

Fills, smoothing, distribution ribbons, mapping views, fragment panels
and heatmap sorting are documented in
[signal-representation.md](signal-representation.md) — all display-layer,
with canonical outputs untouched (smoothing/quantiles export as
separately labeled files; mapping classes are analysis options).

## Reproducibility

Everything above serializes in the config's `style:` section; a
GUI-styled figure reproduces headlessly with `tescope run analysis.yaml`
(and `--style-preset <name>` as a quick override; all presets except `custom`). Configs
without a `style:` section keep working and render exactly as before.

## Example

```yaml
# synthetic/demo labels only
style:
  preset: publication
  palette: colorblind
  fig_width_cm: 17.8
  series:
    control:  {color: "#4477AA"}
    treated:  {color: "#EE6677", line_style: dashed}
    signal_c: {line_width: 2.0, alpha: 0.8}
  show_grid: true
  grid_alpha: 0.25
  legend_position: outside right
  y_min: 0
  title: "Demo figure (synthetic data)"
  heatmap_cmap: magma
```
