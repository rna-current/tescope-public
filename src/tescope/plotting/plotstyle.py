"""Centralized plot-style model: a pure rendering layer.

The style object controls only how figures LOOK — colors, line styles,
typography, legend, display limits, figure size, annotation-lane and
heatmap presentation, export DPI. It is consumed exclusively by the
figure builders (:mod:`tescope.plotting.profiles`) and export helpers;
no analysis, projection, aggregation, normalization or QC code reads it,
so no style option can change a scientific number (pinned by regression
tests).

Presets:
- ``screen`` (default): replicates TEscope's established on-screen look —
  configs without a style section render exactly as before.
- ``publication``: manuscript-oriented — double-column width, restrained
  typography, slightly thinner lines. Generic; no journal is claimed.
- ``publication_minimal``: like ``publication`` with quieter typography
  and thinner lines for dense multi-panel layouts.
- ``presentation``: slide/projector-oriented — wide figure, large type,
  thick lines; Reset additionally selects the filled representation.
- ``high_contrast``: thick lines and enlarged type at screen width for
  low-vision/projector readability.
- ``custom``: same defaults as screen; the label signals user overrides.

Selecting a preset name (config ``style.preset`` / CLI
``--style-preset``) changes only the preset-resolved SIZE values
(typography, line width, figure width); the GUI Reset action /
:meth:`PlotStyleConfig.preset_defaults` additionally applies the
preset's recommended representation defaults. Either way the preset is
pure display — no analysis value can change (pinned by tests).

Display limits (y_min/y_max/x_min/x_max) are VIEW windows only — they
never filter, truncate or transform underlying data, unlike the
analytical filters in the config's ``filters`` section.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

from matplotlib.colors import is_color_like

from .transforms import SMOOTHING_METHODS, TransformError, validate_smoothing

PRESETS = ("screen", "publication", "publication_minimal", "presentation",
           "high_contrast", "custom")

REPRESENTATIONS = ("line", "filled")
BANDS = ("sem", "iqr", "none")
HEATMAP_SORTS = ("cons_start", "genome_length", "consensus_coverage", "total_signal")

LINE_STYLES = ("solid", "dashed", "dotted", "dashdot")

LEGEND_POSITIONS = (
    "best",
    "upper right",
    "upper left",
    "lower right",
    "lower left",
    "outside right",
)

# Curated colormaps for the per-copy heatmap (all matplotlib built-ins,
# perceptually uniform or standard diverging).
HEATMAP_CMAPS = ("viridis", "magma", "cividis", "Greys", "RdBu_r", "coolwarm")

# Curated series palettes. "default" is the established Okabe-Ito cycle
# (already color-vision-conscious); "colorblind" is Paul Tol's bright
# scheme; "grayscale" pairs grays with automatic line-style cycling so
# series stay distinguishable in print.
PALETTES: dict[str, tuple[str, ...]] = {
    "default": (
        "#0072B2", "#D55E00", "#009E73", "#CC79A7",
        "#E69F00", "#56B4E9", "#F0E442", "#000000",
    ),
    "colorblind": (
        "#4477AA", "#EE6677", "#228833", "#CCBB44",
        "#66CCEE", "#AA3377", "#BBBBBB", "#000000",
    ),
    "grayscale": (
        "#000000", "#555555", "#888888", "#AAAAAA",
        "#333333", "#777777", "#999999", "#BBBBBB",
    ),
}

# journal-layout figure-width conveniences (centimetres; approximate
# common single/double column widths — no specific journal is claimed)
FIGURE_WIDTH_PRESETS_CM = {"single-column": 8.5, "double-column": 17.8}

CM_PER_INCH = 2.54


class PlotStyleError(ValueError):
    pass


@dataclass
class SeriesStyle:
    """Per-series visual overrides, keyed by the stable series label
    (group name, or track display name for ungrouped tracks)."""

    color: str = ""  # hex/named color; "" = palette assignment
    line_style: str = ""  # one of LINE_STYLES; "" = preset default
    line_width: float = 0.0  # 0 = preset default
    alpha: float = 1.0
    fill_color: str = ""  # filled representation only; "" = series color

    def validate(self, label: str = "") -> None:
        where = f"series '{label}': " if label else ""
        if self.color and not is_color_like(self.color):
            raise PlotStyleError(f"{where}invalid color {self.color!r}")
        if self.line_style and self.line_style not in LINE_STYLES:
            raise PlotStyleError(
                f"{where}line_style must be one of {LINE_STYLES}, got {self.line_style!r}"
            )
        if self.line_width < 0:
            raise PlotStyleError(f"{where}negative line_width")
        if not (0.0 <= self.alpha <= 1.0):
            raise PlotStyleError(f"{where}alpha must be in [0, 1]")
        if self.fill_color and not is_color_like(self.fill_color):
            raise PlotStyleError(f"{where}invalid fill_color {self.fill_color!r}")

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> SeriesStyle:
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        return cls(**known)


# preset base values; 0.0 / "" fields in PlotStyleConfig mean "use preset"
_PRESET_BASE = {
    "screen": {
        # replicates the pre-style-subsystem rendering exactly
        "fig_width_cm": 7.2 * CM_PER_INCH,
        "title_size": 10.0,
        "label_size": 9.0,
        "tick_size": 8.0,
        "legend_size": 8.0,
        "annotation_size": 6.5,
        "line_width": 1.2,
    },
    "publication": {
        "fig_width_cm": FIGURE_WIDTH_PRESETS_CM["double-column"],
        "title_size": 9.0,
        "label_size": 8.0,
        "tick_size": 7.0,
        "legend_size": 7.0,
        "annotation_size": 6.0,
        "line_width": 1.0,
    },
    "publication_minimal": {
        "fig_width_cm": FIGURE_WIDTH_PRESETS_CM["double-column"],
        "title_size": 8.0,
        "label_size": 7.0,
        "tick_size": 6.5,
        "legend_size": 6.5,
        "annotation_size": 5.5,
        "line_width": 0.9,
    },
    "presentation": {
        "fig_width_cm": 24.0,
        "title_size": 15.0,
        "label_size": 13.0,
        "tick_size": 11.0,
        "legend_size": 12.0,
        "annotation_size": 9.0,
        "line_width": 2.2,
    },
    "high_contrast": {
        "fig_width_cm": 7.2 * CM_PER_INCH,
        "title_size": 12.0,
        "label_size": 11.0,
        "tick_size": 10.0,
        "legend_size": 10.0,
        "annotation_size": 8.0,
        "line_width": 2.4,
    },
}
_PRESET_BASE["custom"] = _PRESET_BASE["screen"]

# Reset-time representation defaults per preset (display-only; applied
# by preset_defaults(), i.e. the GUI Reset action — never silently at
# render time, so existing configs keep their stored appearance).
_PRESET_RESET_OVERLAY: dict[str, dict] = {
    "publication_minimal": {"representation": "line", "show_grid": False},
    "presentation": {"representation": "filled", "fill_alpha": 0.18,
                     "band_alpha": 0.25},
    "high_contrast": {"representation": "line", "band_alpha": 0.3},
}


@dataclass
class PlotStyleConfig:
    """Serializable plot style. All fields optional; empty/zero = preset
    default, so configs without a style section render as before."""

    preset: str = "screen"
    palette: str = "default"
    # figure size (physical centimetres; 0 = automatic/preset)
    fig_width_cm: float = 0.0
    fig_height_cm: float = 0.0  # 0 = automatic per panel count (recommended)
    # typography (points; 0 = preset default)
    title_size: float = 0.0
    label_size: float = 0.0
    tick_size: float = 0.0
    legend_size: float = 0.0
    annotation_size: float = 0.0
    # presentation text ("" = TEscope's default semantics)
    title: str = ""
    x_label: str = ""
    y_label: str = ""  # overlay-mode y-axis label
    # legend (overlay mode)
    show_legend: bool = True
    legend_position: str = "upper right"  # matches the pre-style default look
    # axes / display limits (display-only; never filters data)
    show_grid: bool = False
    grid_alpha: float = 0.3
    y_min: float | None = None
    y_max: float | None = None
    x_min: float | None = None  # consensus-axis view window (bp)
    x_max: float | None = None
    # series
    line_width: float = 0.0  # default for all series; 0 = preset
    series: dict = field(default_factory=dict)  # label -> SeriesStyle dict
    # -- signal representation (display transforms; canonical results are
    # -- NEVER modified — see docs/signal-representation.md) --------------
    representation: str = "line"  # line | filled (line + translucent area)
    fill_alpha: float = 0.15
    # smoothing: DISPLAY-ONLY; unsmoothed profile stays canonical
    smoothing_method: str = "none"  # none|moving_average|gaussian|savgol
    smoothing_window: int = 5  # moving_average / savgol (odd)
    smoothing_sigma: float = 2.0  # gaussian (in bins)
    smoothing_polyorder: int = 2  # savgol
    smoothing_show_raw: bool = True  # raw thin + smoothed strong
    # band: "sem" = uncertainty of the group/track mean (existing SEM);
    # "iqr" = VARIATION ACROSS TE COPIES (25-75%, + 10-90% outer band) —
    # per-copy quantiles, single (ungrouped) continuous tracks only
    band: str = "sem"
    outer_band: bool = True  # add 10-90% ribbon when band == "iqr"
    # feature-region background shading across signal panels
    annotation_shading: bool = False
    annotation_shading_alpha: float = 0.08
    # additional figure compositions (derived display, no new numbers)
    fragment_panels: bool = False  # all vs full-length vs 5'-truncated figure
    # uncertainty band (visibility = existing plot.show_sem; color = series color)
    band_alpha: float = 0.2
    # annotation feature lane (origin stays distinguishable: user lane keeps
    # its hatch and both lanes keep their lane labels regardless of colors)
    annotation_alpha: float = 1.0
    annotation_border_width: float = 0.8
    database_lane_color: str = ""  # "" = default #C7DAE3
    user_lane_color: str = ""  # "" = default #EBDDC3
    # heatmap (display only — never rescales the underlying matrix)
    heatmap_cmap: str = "viridis"
    heatmap_vmin: float | None = None
    heatmap_vmax: float | None = None
    heatmap_colorbar: bool = True
    heatmap_sort: str = "cons_start"  # row ORDER only; values never change
    heatmap_show_profile: bool = False  # aggregate profile panel above heatmap

    # -- validation --------------------------------------------------------

    def validate(self) -> None:
        if self.preset not in PRESETS:
            raise PlotStyleError(f"style.preset must be one of {PRESETS}")
        if self.palette not in PALETTES:
            raise PlotStyleError(f"style.palette must be one of {tuple(PALETTES)}")
        if self.legend_position not in LEGEND_POSITIONS:
            raise PlotStyleError(
                f"style.legend_position must be one of {LEGEND_POSITIONS}"
            )
        if self.heatmap_cmap not in HEATMAP_CMAPS:
            raise PlotStyleError(f"style.heatmap_cmap must be one of {HEATMAP_CMAPS}")
        if self.representation not in REPRESENTATIONS:
            raise PlotStyleError(f"style.representation must be one of {REPRESENTATIONS}")
        if self.band not in BANDS:
            raise PlotStyleError(f"style.band must be one of {BANDS}")
        if self.heatmap_sort not in HEATMAP_SORTS:
            raise PlotStyleError(f"style.heatmap_sort must be one of {HEATMAP_SORTS}")
        if not (0.0 <= self.fill_alpha <= 1.0):
            raise PlotStyleError("style.fill_alpha must be in [0, 1]")
        if not (0.0 <= self.annotation_shading_alpha <= 0.5):
            raise PlotStyleError("style.annotation_shading_alpha must be in [0, 0.5]")
        if self.smoothing_method not in SMOOTHING_METHODS:
            raise PlotStyleError(
                f"style.smoothing_method must be one of {SMOOTHING_METHODS}")
        if self.smoothing_method != "none":
            try:
                validate_smoothing(self.smoothing_method, self.smoothing_window,
                                   self.smoothing_sigma, self.smoothing_polyorder)
            except TransformError as e:
                raise PlotStyleError(str(e)) from e
        for name in ("fig_width_cm", "fig_height_cm", "title_size", "label_size",
                     "tick_size", "legend_size", "annotation_size", "line_width",
                     "annotation_border_width"):
            if getattr(self, name) < 0:
                raise PlotStyleError(f"style.{name} must be >= 0")
        for name in ("grid_alpha", "band_alpha", "annotation_alpha"):
            if not (0.0 <= getattr(self, name) <= 1.0):
                raise PlotStyleError(f"style.{name} must be in [0, 1]")
        for name in ("database_lane_color", "user_lane_color"):
            v = getattr(self, name)
            if v and not is_color_like(v):
                raise PlotStyleError(f"style.{name}: invalid color {v!r}")
        if (self.y_min is not None and self.y_max is not None
                and self.y_max <= self.y_min):
            raise PlotStyleError("style.y_max must be > y_min")
        if (self.x_min is not None and self.x_max is not None
                and self.x_max <= self.x_min):
            raise PlotStyleError("style.x_max must be > x_min")
        for label, sd in (self.series or {}).items():
            if not isinstance(sd, (dict, SeriesStyle)):
                raise PlotStyleError(f"style.series[{label!r}] must be a mapping")
            self.series_style(label).validate(label)

    # -- resolution helpers (renderer-facing) ------------------------------

    def _preset_value(self, name: str):
        return _PRESET_BASE[self.preset if self.preset in _PRESET_BASE else "screen"][name]

    def font(self, which: str) -> float:
        v = getattr(self, f"{which}_size")
        return v if v > 0 else self._preset_value(f"{which}_size")

    def resolved_line_width(self) -> float:
        return self.line_width if self.line_width > 0 else self._preset_value("line_width")

    def fig_width_in(self) -> float:
        w = self.fig_width_cm if self.fig_width_cm > 0 else self._preset_value("fig_width_cm")
        return w / CM_PER_INCH

    def fig_height_in(self, auto_height_in: float) -> float:
        """Explicit height, or the caller's per-panel automatic height."""
        if self.fig_height_cm > 0:
            return self.fig_height_cm / CM_PER_INCH
        # automatic height scales with the width relative to the classic
        # 7.2 in layout so panels keep their aspect at other widths
        return auto_height_in * (self.fig_width_in() / 7.2)

    def series_style(self, label: str) -> SeriesStyle:
        """Stored style for a series label (stale/missing entries are fine)."""
        entry = (self.series or {}).get(label)
        if entry is None:
            return SeriesStyle()
        if isinstance(entry, SeriesStyle):
            return entry
        return SeriesStyle.from_dict(dict(entry))

    def series_color(self, label: str, index: int) -> str:
        """Deterministic color: explicit user choice, else palette by
        display order. A saved explicit color is never reassigned."""
        s = self.series_style(label)
        if s.color:
            return s.color
        pal = PALETTES[self.palette]
        return pal[index % len(pal)]

    def series_line_style(self, label: str, index: int) -> str:
        s = self.series_style(label)
        if s.line_style:
            return s.line_style
        if self.palette == "grayscale":
            # cycle line styles so gray series stay distinguishable in print
            return LINE_STYLES[index % len(LINE_STYLES)]
        return "solid"

    def series_line_width(self, label: str) -> float:
        s = self.series_style(label)
        return s.line_width if s.line_width > 0 else self.resolved_line_width()

    def series_alpha(self, label: str) -> float:
        return self.series_style(label).alpha

    def rc_overrides(self) -> dict:
        """Merged into the base publication rc for this render only."""
        return {
            "axes.titlesize": self.font("title"),
            "axes.labelsize": self.font("label"),
            "xtick.labelsize": self.font("tick"),
            "ytick.labelsize": self.font("tick"),
            "legend.fontsize": self.font("legend"),
            "lines.linewidth": self.resolved_line_width(),
        }

    # -- (de)serialization -------------------------------------------------

    def to_dict(self) -> dict:
        d = asdict(self)
        d["series"] = {
            k: (v.to_dict() if isinstance(v, SeriesStyle) else dict(v))
            for k, v in (self.series or {}).items()
        }
        return d

    @classmethod
    def from_dict(cls, d: dict) -> PlotStyleConfig:
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        obj = cls(**known)
        obj.validate()
        return obj

    @classmethod
    def preset_defaults(cls, preset: str) -> PlotStyleConfig:
        """A fresh style reset to the named preset (the Reset action),
        including the preset's recommended representation defaults."""
        if preset not in PRESETS:
            raise PlotStyleError(f"unknown preset {preset!r}")
        return cls(preset=preset, **_PRESET_RESET_OVERLAY.get(preset, {}))
