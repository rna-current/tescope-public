"""Interactive consensus browser: wheel zoom, navigation, cursor readout.

Display-only interaction over an already-computed AnalysisResult. Nothing
here recomputes or mutates scientific values — zooming changes the shared
x window (all panels are created with ``sharex=True``, so signal,
annotation, support, mapping and sequence panels always show exactly the
same consensus interval), the sequence track re-renders its
zoom-dependent letters, and the cursor readout REPORTS values that
already exist in the result.

Complements (never replaces) the matplotlib navigation toolbar, which
already provides synchronized box zoom, pan, home and back/forward
history on the shared axes.

Coordinates: the axis is TEscope's 0-based continuous consensus axis;
positions REPORTED to the user (cursor, Go to, search) are 1-based
inclusive base numbers (see tescope.core.sequence.coords).
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import QObject, Signal

from ...core.sequence.coords import (
    axis_to_base,
    parse_position_or_range,
)
from ...plotting.sequence_track import refresh_letters

_ZOOM_STEP = 1.25  # wheel factor per notch
_MIN_SPAN = 10.0  # bp; prevents zooming into nothing
_MARKER_COLOR = "#B22222"


class ConsensusNavigator(QObject):
    """Attaches interactive navigation to the current results canvas."""

    cursor_info = Signal(str)  # concise readout for the info label
    view_changed = Signal(float, float)  # visible consensus interval

    def __init__(self, parent=None):
        super().__init__(parent)
        self.canvas = None
        self.result = None
        self._axes: list = []
        self._seq_axes: list = []
        self._cids: list = []
        self._markers: list = []
        self._last_bin = None

    # -- lifecycle ---------------------------------------------------------

    def attach(self, canvas, result) -> None:
        """Bind to a freshly built canvas/figure (detaching any previous)."""
        self.detach()
        self.canvas = canvas
        self.result = result
        fig = canvas.figure
        self._axes = list(fig.axes)
        self._seq_axes = [
            a for a in self._axes if hasattr(a, "_tescope_sequence")
        ]
        self._markers = []
        self._last_bin = None
        self._cids = [
            canvas.mpl_connect("scroll_event", self._on_scroll),
            canvas.mpl_connect("motion_notify_event", self._on_motion),
            canvas.mpl_connect("button_press_event", self._on_click),
            # toolbar pan/box-zoom end in a draw: keep letters in sync
            canvas.mpl_connect("draw_event", self._on_draw),
        ]

    def detach(self) -> None:
        if self.canvas is not None:
            for cid in self._cids:
                try:
                    self.canvas.mpl_disconnect(cid)
                except (AttributeError, RuntimeError):
                    pass
        self.canvas = None
        self.result = None
        self._axes = []
        self._seq_axes = []
        self._cids = []
        self._markers = []

    # -- view control ------------------------------------------------------

    @property
    def cons_length(self) -> int:
        return int(self.result.cons_length) if self.result is not None else 0

    def visible(self) -> tuple[float, float]:
        if not self._axes:
            return (0.0, 0.0)
        lo, hi = self._axes[-1].get_xlim()
        return float(lo), float(hi)

    def set_view(self, x_lo: float, x_hi: float) -> None:
        """Clamped display window on the shared consensus axis."""
        if not self._axes:
            return
        length = float(self.cons_length)
        span = min(max(x_hi - x_lo, _MIN_SPAN), length)
        x_lo = max(0.0, min(x_lo, length - span))
        # sharex: setting one axis moves every aligned panel together
        self._axes[-1].set_xlim(x_lo, x_lo + span)
        self._refresh_sequence()
        self.canvas.draw_idle()
        self.view_changed.emit(x_lo, x_lo + span)

    def full_view(self) -> None:
        if self.result is not None:
            self.set_view(0.0, float(self.cons_length))

    def go_to(self, text: str) -> None:
        """'4127' or '4000-4500' (1-based inclusive). Raises RangeError."""
        b1, b2 = parse_position_or_range(text, self.cons_length)
        if b1 == b2:  # single base: a comfortable window around it
            half = max(_MIN_SPAN * 4, 30.0)
            self.set_view(b1 - 1 - half, b1 + half)
        else:
            self.set_view(b1 - 1, float(b2))

    def zoom_to_bases(self, b1: int, b2: int, margin_frac: float = 0.5) -> None:
        """Zoom to a 1-based inclusive base range with context margin."""
        span = max(b2 - b1 + 1, 1)
        m = max(span * margin_frac, _MIN_SPAN)
        self.set_view(b1 - 1 - m, b2 + m)

    # -- event handlers ----------------------------------------------------

    def _on_scroll(self, event) -> None:
        if event.inaxes is None or self.result is None:
            return
        lo, hi = self.visible()
        span = hi - lo
        factor = 1 / _ZOOM_STEP if event.button == "up" else _ZOOM_STEP
        new_span = span * factor
        # zoom centered near the cursor: keep its data position fixed
        x = event.xdata if event.xdata is not None else (lo + hi) / 2
        frac = (x - lo) / span if span > 0 else 0.5
        self.set_view(x - frac * new_span, x + (1 - frac) * new_span)

    def _on_draw(self, _event) -> None:
        # after toolbar pan/zoom/home, re-evaluate letter visibility once;
        # skip when the window is unchanged (cursor marker redraws etc.)
        for ax in self._seq_axes:
            window = tuple(ax.get_xlim())
            if window == getattr(ax, "_tescope_letter_window", None):
                continue
            ax._tescope_letter_window = window  # noqa: SLF001
            before = len(getattr(ax, "_tescope_letter_artists", []))
            after = refresh_letters(ax)
            if after != before:
                # letters changed: one extra idle draw shows them
                self.canvas.draw_idle()

    def _refresh_sequence(self) -> None:
        for ax in self._seq_axes:
            refresh_letters(ax)

    def _on_motion(self, event) -> None:
        if event.inaxes is None or self.result is None or event.xdata is None:
            return
        res = self.result
        b = axis_to_base(event.xdata, res.cons_length)
        if b is None:
            return
        bin_i = int(
            np.clip(np.searchsorted(res.bin_edges, event.xdata, side="right") - 1,
                    0, len(res.bin_edges) - 2)
        )
        if bin_i == self._last_bin and b == getattr(self, "_last_base", None):
            return  # throttle: only re-render the readout when it changes
        self._last_bin, self._last_base = bin_i, b
        self.cursor_info.emit(self._readout(b, bin_i))

    def _on_click(self, event) -> None:
        """Left click places a position marker line across all panels."""
        if (event.inaxes is None or self.result is None
                or event.xdata is None or event.button != 1
                or self.canvas.toolbar is not None
                and getattr(self.canvas.toolbar, "mode", "")):
            return  # never fight the toolbar's pan/zoom drag modes
        for m in self._markers:
            try:
                m.remove()
            except (ValueError, RuntimeError):
                pass
        self._markers = [
            ax.axvline(event.xdata, color=_MARKER_COLOR, linewidth=0.8,
                       alpha=0.7, zorder=10)
            for ax in self._axes
        ]
        self.canvas.draw_idle()

    # -- readout -----------------------------------------------------------

    def _readout(self, base: int, bin_i: int) -> str:
        """Concise cursor line; only fields that actually exist."""
        res = self.result
        parts = [f"base {base:,} (1-based)"]
        seq = getattr(res, "sequence", None)
        if seq is not None and seq.verified:
            parts.append(f"{seq.sequence.sequence[base - 1]}")
        stat = res.config.aggregation.statistic
        shown = 0
        for t in res.continuous_tracks():
            if t.profile is None or not t.config.visible or shown >= 3:
                continue
            v = t.profile[stat].iloc[bin_i]
            parts.append(
                f"{t.name}: " + (f"{v:.3g}" if np.isfinite(v) else "no data")
            )
            shown += 1
        n_cop = int(res.element_covered.sum(axis=0)[bin_i])
        parts.append(f"copies: {n_cop}")
        feats = getattr(res, "features", None)
        if feats is not None:
            x = base - 1  # feature intervals are 0-based half-open
            names = [
                f"{f.name}" + (" (user)" if f.origin == "user-defined" else "")
                for f in feats.all_features()
                if f.start <= x < f.end
            ]
            if names:
                parts.append("feature: " + ", ".join(names[:2]))
        return "   |   ".join(parts)
