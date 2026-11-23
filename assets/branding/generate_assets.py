#!/usr/bin/env python
"""Regenerate all derived TEscope branding assets from the adopted
vector masters (concept B "profile lens", operator-approved 2026-10).

Masters (authored, vector-native, project-owned — no third-party
artwork, no font files distributed):
  source/tescope_symbol_master.svg  — the lens symbol
  source/tescope_logo_master.svg    — symbol + "TEscope" wordmark

Run from the repo root:  .venv/bin/python assets/branding/generate_assets.py

Derivations (deterministic):
- symbol PNGs: QSvgRenderer (PySide6, offscreen) rasterizes the symbol
  master at each size with a transparent background;
- lockup PNGs: the symbol PNG composited with the wordmark drawn by
  Pillow using the installed Arial Bold (falls back to Liberation Sans
  Bold when Arial is absent; actual face recorded on stdout) —
  "TE" #3C5991, "scope" #0072B2; dark variant swaps the navy text for
  a light tone, symbol unchanged;
- icons: Pillow ICO (256..16) and ICNS (512 base) from the symbol;
  favicon.png/.ico; byte copies into packaging/; runtime copies into
  src/tescope/gui/theme/branding/ (bundled via the PyInstaller spec).
"""
from __future__ import annotations

import os
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
BR = ROOT / "assets/branding"
SYM_SVG = BR / "source/tescope_symbol_master.svg"
ICON_SIZES = (1024, 512, 256, 128, 64, 48, 32, 24, 16)

NAVY = (60, 89, 145, 255)       # TE
ACCENT = (0, 114, 178, 255)     # scope
LIGHT = (220, 228, 242, 255)    # dark-background text tone

FONT_CANDIDATES = (
    "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf",
)


def rasterize_symbol(sizes=ICON_SIZES) -> dict[int, Image.Image]:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QGuiApplication, QImage, QPainter
    from PySide6.QtSvg import QSvgRenderer

    QGuiApplication.instance() or QGuiApplication([])
    renderer = QSvgRenderer(str(SYM_SVG))
    assert renderer.isValid(), SYM_SVG
    out = {}
    for s in sizes:
        img = QImage(s, s, QImage.Format_ARGB32_Premultiplied)
        img.fill(Qt.transparent)
        p = QPainter(img)
        p.setRenderHint(QPainter.Antialiasing)
        renderer.render(p)
        p.end()
        tmp = BR / f".sym_{s}.png"
        img.save(str(tmp))
        out[s] = Image.open(tmp).convert("RGBA")
        tmp.unlink()
    return out


def wordmark_font(px: int):
    for path in FONT_CANDIDATES:
        if Path(path).exists():
            return ImageFont.truetype(path, px), Path(path).name
    raise SystemExit("no Arial/Liberation bold face available")


def build_lockup(symbol: Image.Image, te_color, scope_color,
                 scale: int = 4) -> Image.Image:
    """Mirror the lockup master geometry (viewBox 302x72, symbol at
    (4,4) size 64, text baseline y=47 size 34) at `scale`x."""
    W, H = 302 * scale, 72 * scale
    img = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    sym = symbol.resize((64 * scale, 64 * scale), Image.LANCZOS)
    img.paste(sym, (4 * scale, 4 * scale), sym)
    font, _face = wordmark_font(34 * scale)
    d = ImageDraw.Draw(img)
    x = 84 * scale
    baseline = 47 * scale
    asc, _desc = font.getmetrics()
    y = baseline - asc
    d.text((x, y), "TE", font=font, fill=te_color)
    x += d.textlength("TE", font=font) + scale  # ~0.5px letter-spacing*2
    d.text((x, y), "scope", font=font, fill=scope_color)
    end = x + d.textlength("scope", font=font)
    return img.crop((0, 0, min(W, int(end) + 6 * scale), H))


def main() -> None:
    logo_d, sym_d, icons_d = BR / "logo", BR / "symbol", BR / "icons"
    gui_d = ROOT / "src/tescope/gui/theme/branding"
    for d in (logo_d, sym_d, icons_d, gui_d):
        d.mkdir(parents=True, exist_ok=True)

    syms = rasterize_symbol()
    base = syms[1024]
    base.resize((303, 303), Image.LANCZOS).save(sym_d / "tescope_symbol.png")
    for s in (512, 256, 128, 64, 48, 32, 16):
        syms[s].save(sym_d / f"tescope_symbol_{s}.png")

    lock2x = build_lockup(base, NAVY, ACCENT, scale=4)
    main2x = lock2x.resize((lock2x.width // 2, lock2x.height // 2),
                           Image.LANCZOS)
    main1x = lock2x.resize((lock2x.width // 4, lock2x.height // 4),
                           Image.LANCZOS)
    main2x.save(logo_d / "tescope_logo_main@2x.png")
    main1x.save(logo_d / "tescope_logo_main.png")
    main1x.save(logo_d / "tescope_logo_light.png")
    dark2x = build_lockup(base, LIGHT, (86, 180, 233, 255), scale=4)
    dark2x.resize((dark2x.width // 2, dark2x.height // 2),
                  Image.LANCZOS).save(logo_d / "tescope_logo_dark@2x.png")
    dark2x.resize((dark2x.width // 4, dark2x.height // 4),
                  Image.LANCZOS).save(logo_d / "tescope_logo_dark.png")

    syms[256].save(icons_d / "tescope.ico",
                   sizes=[(s, s) for s in (256, 128, 64, 48, 32, 24, 16)])
    syms[512].save(icons_d / "tescope.icns")
    syms[48].save(icons_d / "favicon.png")
    syms[48].save(icons_d / "favicon.ico",
                  sizes=[(16, 16), (32, 32), (48, 48)])
    (ROOT / "packaging/tescope.ico").write_bytes(
        (icons_d / "tescope.ico").read_bytes())
    (ROOT / "packaging/tescope.icns").write_bytes(
        (icons_d / "tescope.icns").read_bytes())

    for s in (16, 24, 32, 48, 64, 128, 256, 512):
        syms[s].save(gui_d / f"tescope_symbol_{s}.png")
    main1x.save(gui_d / "tescope_logo_main.png")

    print("regenerated profile-lens branding;",
          f"lockup {main1x.size} (1x) / {main2x.size} (2x);",
          "wordmark face:", wordmark_font(10)[1])


if __name__ == "__main__":
    main()
