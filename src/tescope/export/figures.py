"""Figure export (PNG/SVG/PDF) with reproducibility metadata.

Every exported figure gets a JSON sidecar recording the TEscope version,
the full analysis config, and annotation provenance — enough to
regenerate the figure from scratch.
"""

from __future__ import annotations

import json
from pathlib import Path

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from ..core.analysis import AnalysisResult
from ..plotting.style import publication_style

VALID_FIGURE_FORMATS = ("png", "svg", "pdf", "tiff")

# export conveniences (no journal-specific requirements hard-coded):
# vector formats ignore DPI by nature; raster presets choose one.
EXPORT_PRESETS = {
    "publication-vector": {"formats": ["pdf", "svg"], "dpi": None},
    "publication-raster": {"formats": ["png", "tiff"], "dpi": 600},
    "screen": {"formats": ["png"], "dpi": 150},
}


def save_figure(
    fig: Figure,
    stem: str | Path,
    formats: list[str],
    result: AnalysisResult | None = None,
    dpi: int | None = None,
) -> list[Path]:
    """Save ``fig`` as ``<stem>.<fmt>`` for each format; returns paths written."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    # Headless figures need a rendering canvas; a figure already embedded in
    # a GUI canvas (Qt) must keep it — replacing it would detach the figure
    # from the widget displaying it.
    if type(fig.canvas).__name__ in ("FigureCanvasBase", "NoneType"):
        FigureCanvasAgg(fig)
    written: list[Path] = []
    if dpi is None:
        dpi = result.config.plot.dpi if result else 300
    # An outside-positioned legend needs a tight bounding box to avoid
    # clipping; the figure builder marks that case. Vector formats keep
    # true vector output (no silent rasterization); dpi applies to raster.
    bbox = "tight" if getattr(fig, "_tescope_tight_bbox", False) else None
    # publication_style() at SAVE time guarantees the backend export
    # settings (Arial-first font.family metadata, svg.fonttype "none",
    # pdf.fonttype 42) apply even when a caller saves outside the
    # rendering context — text geometry was already fixed at build time.
    with publication_style():
        for fmt in formats:
            if fmt not in VALID_FIGURE_FORMATS:
                raise ValueError(f"unsupported figure format {fmt!r}")
            p = stem.with_suffix(f".{fmt}")
            fig.savefig(p, format=fmt, dpi=dpi, bbox_inches=bbox)
            written.append(p)
    if result is not None:
        meta = {
            "tescope_version": _version(),
            "figure": stem.name,
            "config": result.config.to_dict(),
            "annotation_provenance": result.provenance,
            "n_annotated": result.n_annotated,
            "n_selected": result.n_selected,
        }
        seq = getattr(result, "sequence", None)
        if seq is not None:
            meta["sequence_provenance"] = {
                "status": seq.status,
                "detail": seq.detail,
                **(seq.provenance.to_dict() if seq.provenance else {}),
            }
        mp = stem.with_suffix(".figure.json")
        mp.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        written.append(mp)
    return written


def render_result_figures(result: AnalysisResult, out_dir: str | Path) -> list[Path]:
    """Render and save the standard figures for an analysis."""
    from ..plotting.profiles import build_figure, build_heatmap

    out = Path(out_dir)
    formats = result.config.output.figure_formats
    written: list[Path] = []
    with publication_style():
        fig = build_figure(result, mode="stacked")
        written += save_figure(fig, out / "profile_stacked", formats, result)
        if len(list(result.continuous_tracks())) + len(result.group_profiles) > 1:
            fig = build_figure(result, mode="overlay")
            written += save_figure(fig, out / "profile_overlay", formats, result)
        style = getattr(result.config, "style", None)
        if style is not None and style.fragment_panels and \
                "completeness" in result.copies.columns:
            fig = build_figure(result, mode="fragments")
            written += save_figure(fig, out / "profile_fragments", formats, result)
        for t in result.continuous_tracks():
            if t.matrix is not None and t.matrix.n_copies >= 2:
                hm = build_heatmap(result, t)
                written += save_figure(
                    hm, out / f"heatmap_{t.name.replace('/', '_')}", formats, result
                )
    return written


def _version() -> str:
    from tescope import __version__

    return __version__
