"""Consensus sequence track: zoom-dependent, coordinate-exact rendering.

Drawn ONLY for a VERIFIED sequence (resolution guarantees model
compatibility). Three representations depending on the visible span and
available pixel width — never thousands of unreadable letters:

- full/wide view: a compact per-base color band (single imshow artist —
  efficient at any length) acting as a presence/availability strip;
- intermediate zoom: the same band (colors resolve into visible blocks);
- base-level zoom (few enough visible bases for readable glyphs): actual
  A/C/G/T/N letters aligned to their positions, on top of the band.

Colors are restrained and never the only carrier of identity: at
sufficient zoom the letter itself is visible (grayscale-safe); the band
alone is only ever an overview. Letter artists are created for the
visible range only and pruned on every update (performance: a handful of
text artists, not cons_length of them).

Base i (1-based) occupies axis interval [i-1, i); letters sit at i-0.5.
"""

from __future__ import annotations

import numpy as np
from matplotlib.colors import ListedColormap

# muted per-base colors (identity is carried by the letter at zoom;
# the band is an overview). Index order: A C G T U N other-IUPAC.
_BASE_ORDER = "ACGTUN"
_BASE_COLORS = [
    "#7FBF7F",  # A - muted green
    "#7F9FD9",  # C - muted blue
    "#D9B87F",  # G - muted amber
    "#D98A8A",  # T - muted red
    "#D98A8A",  # U - as T
    "#C8C8C8",  # N - gray
    "#BDBDBD",  # other IUPAC ambiguity codes
]
_CMAP = ListedColormap(_BASE_COLORS)

# a letter needs at least this many horizontal pixels to be readable
_MIN_PX_PER_BASE = 7.0
# hard cap on letter artists regardless of pixel math
_MAX_LETTERS = 400


def base_indices(seq: str) -> np.ndarray:
    """Sequence string -> color indices for the band imshow."""
    lut = np.full(256, len(_BASE_ORDER), dtype=np.uint8)
    for i, b in enumerate(_BASE_ORDER):
        lut[ord(b)] = i
    return lut[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def draw_sequence_band(ax, seq: str, cons_length: int) -> None:
    """The static per-base color band across the full consensus."""
    idx = base_indices(seq)
    ax.imshow(
        idx[np.newaxis, :],
        aspect="auto",
        interpolation="nearest",
        extent=(0, cons_length, 0, 1),
        cmap=_CMAP,
        vmin=0,
        vmax=len(_BASE_COLORS) - 1,
        zorder=1,
    )
    ax.set_yticks([])
    ax.set_ylim(0, 1)
    ax.set_ylabel("consensus\nsequence", rotation=0, ha="right", va="center")
    for side in ("left",):
        ax.spines[side].set_visible(False)
    # marks this axis for the interactive navigator and letter refresh
    ax._tescope_sequence = seq  # noqa: SLF001 - own attribute


def _axis_width_px(ax) -> float:
    try:
        bbox = ax.get_window_extent()
        return float(bbox.width)
    except (RuntimeError, AttributeError):
        # no renderer yet: estimate from figure size and dpi
        fig = ax.figure
        return fig.get_figwidth() * fig.dpi * 0.8


def letters_visible(ax, x_lo: float, x_hi: float) -> bool:
    """Whether base letters would be readable for this window."""
    span = max(x_hi - x_lo, 1e-9)
    if span > _MAX_LETTERS:
        return False
    return _axis_width_px(ax) / span >= _MIN_PX_PER_BASE


def refresh_letters(ax, x_lo: float | None = None,
                    x_hi: float | None = None) -> int:
    """(Re)draw base letters for the visible window; returns letters drawn.

    Removes previous letter artists first. Above the readability
    threshold no letters are drawn (the band remains). Shared by the
    interactive navigator (on every zoom/pan) and the headless export
    path (once, for the configured view window).
    """
    seq = getattr(ax, "_tescope_sequence", None)
    if seq is None:
        return 0
    for t in getattr(ax, "_tescope_letter_artists", []):
        t.remove()
    ax._tescope_letter_artists = []  # noqa: SLF001

    if x_lo is None or x_hi is None:
        x_lo, x_hi = ax.get_xlim()
    if not letters_visible(ax, x_lo, x_hi):
        return 0
    lo = max(int(np.floor(x_lo)), 0)
    hi = min(int(np.ceil(x_hi)), len(seq))
    if hi <= lo:
        return 0
    span = max(x_hi - x_lo, 1e-9)
    px_per_base = _axis_width_px(ax) / span
    fontsize = float(np.clip(px_per_base * 0.9, 6.0, 11.0))
    artists = []
    for b0 in range(lo, hi):  # b0 = 0-based index; letter of base b0+1
        artists.append(
            ax.text(
                b0 + 0.5, 0.5, seq[b0],
                ha="center", va="center",
                fontsize=fontsize,
                family=["DejaVu Sans Mono", "monospace"],  # mpl-bundled first
                color="#1A1A1A", zorder=3, clip_on=True,
            )
        )
    ax._tescope_letter_artists = artists  # noqa: SLF001
    return len(artists)
