"""Publication plot style. Applied explicitly, never globally on import.

Typography policy (2026-09): **Arial-first**. Every piece of figure
text (titles, axis/tick labels, legends, annotations, support-panel
and feature-lane labels) uses the declared family chain
``Arial -> Liberation Sans -> DejaVu Sans -> installed CJK fallbacks
-> sans-serif``. Hosts without native Arial render with the first
available family (Liberation Sans is metrically Arial-compatible);
the declared preference stays Arial everywhere, including the
font-family metadata of SVG exports (svg.fonttype "none"). PDF embeds
the RESOLVED face as TrueType (fonttype 42) — i.e. real Arial where
installed, otherwise the first available fallback; this backend
limitation is documented rather than worked around by rasterizing.
"""

from __future__ import annotations

import logging
import re
from contextlib import contextmanager

import matplotlib as mpl

# Declared Latin preference chain — kept in font.family UNCONDITIONALLY
# so exports always carry the intended Arial-first declaration, even on
# hosts where Arial is absent. NOTE: the list must be font.family
# itself — the "sans-serif" alias + font.sans-serif list does NOT build
# a fallback chain (verified on mpl 3.11).
PREFERRED_FAMILIES = ("Arial", "Liberation Sans", "DejaVu Sans")

# Glyph-level CJK/etc. fallback (matplotlib >= 3.6): user-provided
# labels (track names, repeat names) containing non-Latin glyphs fall
# back to whichever of these families the host provides instead of
# rendering as boxes.
_FALLBACK_CANDIDATES = (
    "Noto Sans KR", "Noto Sans CJK KR", "Noto Sans CJK JP",
    "Noto Sans CJK SC", "Malgun Gothic", "Yu Gothic",
    "Hiragino Sans", "Microsoft YaHei", "PingFang SC",
    "Arial Unicode MS",
)


class _MissingPreferredFontFilter(logging.Filter):
    """Silence matplotlib's per-draw findfont warnings ONLY for our
    deliberately-declared preferred families (Arial may legitimately be
    absent, e.g. on Linux); every other font warning passes through."""

    _families = "|".join(re.escape(f) for f in PREFERRED_FAMILIES)
    _pat = re.compile(
        f"findfont: Font family (['\"])({_families})\\1 not found")

    def filter(self, record: logging.LogRecord) -> bool:  # pragma: no cover
        return not self._pat.search(record.getMessage())


logging.getLogger("matplotlib.font_manager").addFilter(
    _MissingPreferredFontFilter())


def _installed_fallbacks() -> list[str]:
    """CJK candidates actually present on this host — matplotlib logs a
    findfont warning for every absent family on every draw, so this
    (long) optional tail is filtered once here."""
    try:
        from matplotlib.font_manager import fontManager

        names = {f.name for f in fontManager.ttflist}
    except Exception:  # pragma: no cover - degenerate env
        return []
    return [f for f in _FALLBACK_CANDIDATES if f in names]


def font_family_chain() -> list[str]:
    """The full declared family chain used for every figure."""
    return [*PREFERRED_FAMILIES, *_installed_fallbacks(), "sans-serif"]


STYLE = {
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "font.family": font_family_chain(),
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "lines.linewidth": 1.2,
    "svg.fonttype": "none",  # editable text in SVG
    "pdf.fonttype": 42,  # TrueType in PDF (editable in Illustrator)
}

# Colorblind-friendly default cycle (Okabe-Ito)
COLOR_CYCLE = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#F0E442", "#000000",
]


@contextmanager
def publication_style():
    with mpl.rc_context(STYLE):
        yield


def color_for(index: int, override: str = "") -> str:
    return override or COLOR_CYCLE[index % len(COLOR_CYCLE)]
