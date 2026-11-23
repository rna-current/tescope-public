"""TEscope icons: approved branding assets + original SVG fallbacks.

The application icon and logo come from the operator-approved TEscope
branding artwork (assets/branding/, provenance in
assets/branding/BRANDING.md); raster copies for runtime use live in
``theme/branding/`` and are bundled into frozen builds via the
PyInstaller spec's ``datas`` entry. The small track-kind icons remain
the original in-memory SVG set (MIT, project-owned). If the branding
files are missing, the app icon falls back to the previous original
SVG mark so startup never fails on a missing asset.
"""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import QByteArray, Qt
from PySide6.QtGui import QIcon, QImage, QPainter, QPixmap
from PySide6.QtSvg import QSvgRenderer

# -- approved branding assets ----------------------------------------------

_BRANDING_SIZES = (16, 24, 32, 48, 64, 128, 256, 512)


def _branding_dir() -> Path:
    """Locate the packaged branding assets in source and frozen runs."""
    if getattr(sys, "frozen", False):  # PyInstaller bundle
        return (Path(getattr(sys, "_MEIPASS", "."))
                / "tescope/gui/theme/branding")
    return Path(__file__).resolve().parent / "branding"


def branding_logo_pixmap() -> QPixmap:
    """The approved main wordmark logo (empty pixmap if unavailable)."""
    path = _branding_dir() / "tescope_logo_main.png"
    return QPixmap(str(path)) if path.exists() else QPixmap()


# -- app icon ---------------------------------------------------------------

_APP_SVG = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64">
  <rect x="1" y="1" width="62" height="62" rx="13" fill="{bg}"/>
  <!-- profile trace -->
  <path d="M10 26 C 16 26, 18 12, 24 12 S 32 26, 38 26 S 48 18, 54 20"
        fill="none" stroke="{trace}" stroke-width="3.4" stroke-linecap="round"/>
  <!-- consensus + aligned fragments -->
  <rect x="10" y="34" width="44" height="4.5" rx="2.25" fill="{cons}"/>
  <rect x="18" y="43" width="36" height="4.5" rx="2.25" fill="{frag}"/>
  <rect x="34" y="52" width="20" height="4.5" rx="2.25" fill="{frag}"/>
</svg>"""


def app_icon(dark: bool = False) -> QIcon:
    """Approved TEscope symbol (lens/DNA mark); ``dark`` kept for API
    compatibility — the branded symbol works on both themes. Falls back
    to the original SVG mark if the branding assets are absent."""
    bdir = _branding_dir()
    icon = QIcon()
    for s in _BRANDING_SIZES:
        p = bdir / f"tescope_symbol_{s}.png"
        if p.exists():
            icon.addPixmap(QPixmap(str(p)))
    if not icon.isNull():
        return icon
    svg = _APP_SVG.format(
        bg="#166E8F" if not dark else "#134F66",
        trace="#FFFFFF",
        cons="#FFFFFF",
        frag="#9FD1E2",
    )
    return _svg_icon(svg, sizes=(16, 24, 32, 48, 64, 128, 256))


# -- track-kind icons -------------------------------------------------------

_KIND_SVGS = {
    "continuous": """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16">
  <path d="M1.5 11 C 3.5 11, 4 4.5, 6 4.5 S 8.5 11, 10.5 11 S 13 7, 14.5 7.6"
        fill="none" stroke="{c}" stroke-width="1.8" stroke-linecap="round"/>
</svg>""",
    "interval": """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16">
  <rect x="1.5" y="4" width="7" height="3.2" rx="1.1" fill="{c}"/>
  <rect x="6.5" y="9" width="8" height="3.2" rx="1.1" fill="{c}"/>
</svg>""",
    "point": """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16">
  <circle cx="3.5" cy="10" r="1.9" fill="{c}"/>
  <circle cx="8" cy="5.5" r="1.9" fill="{c}"/>
  <circle cx="12.5" cy="9" r="1.9" fill="{c}"/>
</svg>""",
}

_ICON_CACHE: dict[tuple[str, str], QIcon] = {}


def kind_icon(kind: str, color: str = "#6A7076") -> QIcon:
    """Small badge icon for a track kind ('continuous'/'interval'/'point')."""
    key = (kind, color)
    if key not in _ICON_CACHE:
        svg = _KIND_SVGS.get(kind)
        if svg is None:
            return QIcon()
        _ICON_CACHE[key] = _svg_icon(svg.format(c=color), sizes=(16, 24, 32))
    return _ICON_CACHE[key]


# -- helpers ----------------------------------------------------------------


def _svg_icon(svg: str, sizes: tuple[int, ...]) -> QIcon:
    """Render the SVG crisply at each requested size."""
    icon = QIcon()
    renderer = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    if not renderer.isValid():
        return icon
    for s in sizes:
        img = QImage(s, s, QImage.Format_ARGB32_Premultiplied)
        img.fill(Qt.transparent)
        painter = QPainter(img)
        painter.setRenderHint(QPainter.Antialiasing)
        renderer.render(painter)
        painter.end()
        icon.addPixmap(QPixmap.fromImage(img))
    return icon
