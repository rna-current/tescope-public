#!/usr/bin/env python
"""Generate the TEscope User Guide cover (SVG source of truth).

Inlines the vector logo master (assets/branding/source/
tescope_logo_master.svg) into an A4 portrait cover — self-contained
SVG, no external references. PDF/PNG renders are produced by the
commands documented in ../README.md (headless Chrome).

Run from repo root:
  .venv/bin/python docs/user-guide/assets/make_cover.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
MASTER = ROOT / "assets/branding/source/tescope_logo_master.svg"
OUT = Path(__file__).resolve().parent

W, H = 794, 1123  # A4 portrait @ 96 dpi
NAVY = "#3C5991"
GREY = "#5A6472"


def main() -> None:
    svg = MASTER.read_text()
    m = re.search(r"<defs>.*?</defs>", svg, re.S)
    defs = m.group(0) if m else ""
    inner = svg.split("</defs>")[1] if m else re.search(
        r"<svg[^>]*>(.*)</svg>", svg, re.S).group(1)
    body = inner.rsplit("</svg>")[0]
    # logo viewBox is (0,0) 302x72 -> place at width 360
    scale = 360 / 302
    tx = (W - 360) / 2
    ty = 165
    commit = subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
        capture_output=True, text=True).stdout.strip() or "unknown"
    sys.path.insert(0, str(ROOT / "src"))
    from tescope import __version__ as version

    cover = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">
{defs}
  <rect width="{W}" height="{H}" fill="#FFFFFF"/>
  <g transform="translate({tx:.1f} {ty:.1f}) scale({scale:.5f})">
{body}
  </g>
  <g font-family="'Liberation Sans', Helvetica, Arial, sans-serif"
     text-anchor="middle">
    <text x="{W / 2}" y="470" font-size="40" font-weight="bold"
          fill="{NAVY}">TEscope User Guide</text>
    <text x="{W / 2}" y="512" font-size="17" fill="{GREY}">
      Fragment-aware consensus-coordinate profiling of
    </text>
    <text x="{W / 2}" y="536" font-size="17" fill="{GREY}">
      transposable-element signals
    </text>
    <line x1="{W / 2 - 120}" y1="580" x2="{W / 2 + 120}" y2="580"
          stroke="#C9D2DE" stroke-width="1"/>
    <text x="{W / 2}" y="625" font-size="14" fill="{GREY}">
      For installation, usage, and reproducibility guidance
    </text>
    <text x="{W / 2}" y="1020" font-size="12" fill="{GREY}">
      TEscope {version}
    </text>
    <text x="{W / 2}" y="1042" font-size="12" fill="{GREY}">
      https://github.com/rna-current/tescope-public
    </text>
  </g>
</svg>
"""
    (OUT / "cover.svg").write_text(cover)
    (OUT / "cover_print.html").write_text(
        "<!doctype html><html><head><style>"
        "@page{size:A4;margin:0}html,body{margin:0;padding:0}"
        "img{width:210mm;height:297mm;display:block}"
        "</style></head><body>"
        '<img src="cover.svg"></body></html>\n')
    print("wrote cover.svg + cover_print.html (logo inlined from the "
          f"vector master; commit {commit})")


if __name__ == "__main__":
    main()
