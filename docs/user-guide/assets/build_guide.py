#!/usr/bin/env python3
"""Build the single-file TEscope User Guide (HTML + PDF).

Concatenates cover + README.md + reference.md + concepts.md into one
self-contained HTML file (all images inlined as data URIs), then
prints it to PDF with headless Chrome.

The ``markdown`` package is NOT a project dependency; point
``GUIDE_PYDEPS`` at a directory containing it (e.g. installed with
``uv pip install --target <dir> markdown``) or have it importable.

Usage:  python docs/user-guide/assets/build_guide.py
Output: docs/user-guide/TEscope_User_Guide.html / .pdf
"""

from __future__ import annotations

import base64
import os
import re
import subprocess
import sys
from pathlib import Path

_pydeps = os.environ.get("GUIDE_PYDEPS")
if _pydeps:
    sys.path.insert(0, _pydeps)
import markdown  # noqa: E402

GUIDE = Path(__file__).resolve().parents[1]
REPO = GUIDE.parents[2]
sys.path.insert(0, str(REPO / "src"))
from tescope import __version__  # noqa: E402

PARTS = [
    ("quickstart", GUIDE / "README.md"),
    ("reference", GUIDE / "reference.md"),
    ("concepts", GUIDE / "concepts.md"),
]

MIME = {".png": "image/png", ".svg": "image/svg+xml", ".jpg": "image/jpeg"}


def data_uri(path: Path) -> str:
    mime = MIME[path.suffix.lower()]
    return f"data:{mime};base64,{base64.b64encode(path.read_bytes()).decode()}"


def md_to_html(md_path: Path) -> str:
    text = md_path.read_text(encoding="utf-8")
    html = markdown.markdown(text, extensions=["tables", "fenced_code", "sane_lists"])
    base = md_path.parent

    def inline_img(m: re.Match[str]) -> str:
        src = m.group(1)
        p = (base / src).resolve()
        if p.is_file() and p.suffix.lower() in MIME:
            return m.group(0).replace(src, data_uri(p))
        raise SystemExit(f"missing image: {src} (from {md_path.name})")

    html = re.sub(r'<img[^>]*?src="([^":]+)"', inline_img, html)

    # Links between the three bundled parts -> in-document anchors.
    html = re.sub(r'href="reference\.md(#[^"]*)?"', 'href="#reference"', html)
    html = re.sub(r'href="concepts\.md(#[^"]*)?"', 'href="#concepts"', html)

    # Links to repo docs/assets can't resolve in a standalone file:
    # keep the text, show the repository path.
    def repo_link(m: re.Match[str]) -> str:
        return f"<code>docs/{Path(m.group(1)).name}</code>"

    html = re.sub(r'<a href="\.\./([^"]+\.md)">(.+?)</a>', repo_link, html)
    html = re.sub(
        r'<a href="(?:assets/[^"]+|\.\./\.\./[^"]+)">(.+?)</a>',
        r"\1",
        html,
    )
    return html


CSS = """
body { font-family: Arial, 'Liberation Sans', 'DejaVu Sans', sans-serif;
       max-width: 52rem; margin: 0 auto; padding: 0 1.2rem 3rem;
       color: #1a1a1a; line-height: 1.45; font-size: 10.5pt; }
h1 { font-size: 17pt; border-bottom: 2px solid #0b5d4f;
     padding-bottom: .2em; margin-top: 2em; }
h2 { font-size: 13pt; color: #0b5d4f; margin-top: 1.6em; }
h3 { font-size: 11pt; }
table { border-collapse: collapse; width: 100%; margin: .7em 0;
        font-size: 9.5pt; }
th, td { border: 1px solid #bbb; padding: .3em .5em;
         text-align: left; vertical-align: top; }
th { background: #eef4f2; }
code { font-family: 'DejaVu Sans Mono', Consolas, monospace;
       font-size: 90%; background: #f4f4f4; padding: 0 .2em; }
pre { background: #f4f4f4; padding: .6em; overflow-x: auto;
      font-size: 9pt; }
pre code { background: none; padding: 0; }
img { max-width: 100%; border: 1px solid #ddd; }
.cover { text-align: center; page-break-after: always; padding-top: 8em; }
.cover img { border: none; }
.part { page-break-before: always; }
tr, pre, img { page-break-inside: avoid; }
@media print { body { max-width: none; } }
"""


def main() -> None:
    cover = GUIDE / "assets" / "cover.svg"
    sections = []
    sections.append(
        '<div class="cover">'
        f'<img src="{data_uri(cover)}" alt="TEscope User Guide cover" '
        'style="width:70%"></div>'
    )
    for anchor, path in PARTS:
        sections.append(f'<div class="part" id="{anchor}">{md_to_html(path)}</div>')

    html = (
        "<!DOCTYPE html>\n<html lang='en'><head><meta charset='utf-8'>"
        f"<title>TEscope {__version__} User Guide</title>"
        f"<style>{CSS}</style></head><body>" + "\n".join(sections) + "</body></html>"
    )
    out_html = GUIDE / "TEscope_User_Guide.html"
    out_html.write_text(html, encoding="utf-8")
    print(f"wrote {out_html} ({out_html.stat().st_size / 1e6:.1f} MB)")

    out_pdf = GUIDE / "TEscope_User_Guide.pdf"
    chrome = os.environ.get("GUIDE_CHROME", "google-chrome")
    try:
        subprocess.run(
            [
                chrome,
                "--headless",
                "--disable-gpu",
                f"--print-to-pdf={out_pdf}",
                "--no-pdf-header-footer",
                f"file://{out_html}",
            ],
            check=True,
            capture_output=True,
            timeout=300,
        )
        print(f"wrote {out_pdf} ({out_pdf.stat().st_size / 1e6:.1f} MB)")
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"PDF step skipped/failed ({exc}); HTML build is complete.")


if __name__ == "__main__":
    main()
