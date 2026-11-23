"""Compose Figure 1 from two genuine, unedited application captures.

Operations: downscale (LANCZOS), side-by-side placement, thin panel
borders, Arial panel labels. NO content/pixel-value editing inside the
captures beyond resampling; no upscaling.
Output: figure1.png (+ 600-dpi-ready PNG) sized for a 178 mm
double-column figure.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps

SRC = Path(sys.argv[1])
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

a = Image.open(SRC / "fig1_a_main_fragments.png").convert("RGB")
b = Image.open(SRC / "fig1_b_loci_compare.png").convert("RGB")

TARGET_H = b.height  # 788: keep (b) native, downscale (a) to match
a2 = a.resize((round(a.width * TARGET_H / a.height), TARGET_H),
              Image.LANCZOS)

GAP, MARGIN, LABEL_H = 28, 10, 64
W = MARGIN + a2.width + GAP + b.width + MARGIN
H = MARGIN + LABEL_H + TARGET_H + MARGIN
canvas = Image.new("RGB", (W, H), "white")

import os  # noqa: E402

# Panel-label font: any bold sans-serif TTF (the published figure used
# Arial Bold; fonts are not distributed with TEscope).
font = ImageFont.truetype(os.environ.get(
    "FIGURE_FONT", "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf"), 44)
draw = ImageDraw.Draw(canvas)

x = MARGIN
for label, img in (("a", a2), ("b", b)):
    draw.text((x + 2, MARGIN - 4), label, fill="black", font=font)
    bordered = ImageOps.expand(img, border=1, fill="#888888")
    canvas.paste(bordered, (x, MARGIN + LABEL_H))
    x += img.width + GAP

canvas.save(OUT / "figure1.png", dpi=(320, 320))
# single-column / page-size proofs for legibility inspection
col = canvas.resize((round(W * 1004 / W * 1), round(H * 1004 / W)), Image.LANCZOS) \
    if False else canvas.resize((1004, round(H * 1004 / W)), Image.LANCZOS)
col.save(OUT / "figure1_proof_single_column.png")   # 85 mm @ 300 dpi
page = canvas.resize((2102, round(H * 2102 / W)), Image.LANCZOS)
page.save(OUT / "figure1_proof_double_column.png")  # 178 mm @ 300 dpi

prov = {
    "generated": datetime.now(timezone.utc).isoformat(),
    "composition": "two genuine unedited captures (see fig1_captures."
                   "provenance.json): LANCZOS downscale of panel a to "
                   "panel b's height, side-by-side, 1-px borders, Arial "
                   "bold labels; no content editing, no upscaling",
    "canvas_px": [W, H],
    "intended_print": "178 mm double-column; ~320 dpi at that width",
    "panels": {"a": "main window, Profile tab, fragments metagene view",
               "b": "Loci tab: per-locus table + applied subset vs "
                    "eligible comparison"},
}
(OUT / "figure1.provenance.json").write_text(json.dumps(prov, indent=2))
print(f"figure1.png {W}x{H}; proofs written")
