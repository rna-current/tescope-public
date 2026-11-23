"""Base-coordinate conventions — the one place they are converted.

TEscope's internal consensus axis is 0-based half-open (bin edges,
feature intervals, sequence string indices). USER-FACING base positions
(cursor readout, Go to, search hits, selections, FASTA headers) are
1-BASED INCLUSIVE, the convention of sequence tools. The conversion is
implemented and tested here explicitly; nothing else converts.

Base b (1-based) occupies the axis interval [b-1, b).
"""

from __future__ import annotations

import re


class RangeError(ValueError):
    pass


def axis_to_base(x: float, cons_length: int) -> int | None:
    """1-based base number at continuous axis position ``x``.

    None outside [0, cons_length) — positions past the end are not a base.
    """
    if x < 0 or x >= cons_length:
        return None
    return int(x) + 1


def base_to_axis_span(base_start: int, base_end: int) -> tuple[float, float]:
    """1-based inclusive base range -> 0-based half-open axis interval."""
    return float(base_start - 1), float(base_end)


def base_range_to_slice(base_start: int, base_end: int) -> slice:
    """1-based inclusive base range -> Python slice of the sequence string."""
    return slice(base_start - 1, base_end)


def visible_base_range(x_lo: float, x_hi: float, cons_length: int) -> tuple[int, int]:
    """1-based inclusive base range of bases fully or partly visible in an
    axis window, clamped to the sequence."""
    lo = max(int(x_lo), 0) + 1
    hi = min(int(-(-x_hi // 1)), cons_length)  # ceil, clamped
    if hi < lo:
        hi = lo
    return lo, int(hi)


_RANGE_RE = re.compile(r"^\s*(\d[\d,]*)\s*(?:[-:]\s*(\d[\d,]*)\s*)?$")


def parse_position_or_range(text: str, cons_length: int) -> tuple[int, int]:
    """Parse '4127' or '4000-4500' (1-based inclusive), validating bounds.

    A single position returns (pos, pos). Raises RangeError with the
    valid bounds on any invalid input.
    """
    m = _RANGE_RE.match(text or "")
    if not m:
        raise RangeError(
            "enter a base position like 4127 or a range like 4000-4500 "
            "(1-based)"
        )
    b1 = int(m.group(1).replace(",", ""))
    b2 = int(m.group(2).replace(",", "")) if m.group(2) else b1
    if b2 < b1:
        b1, b2 = b2, b1
    if b1 < 1 or b2 > cons_length:
        raise RangeError(
            f"positions must be within 1-{cons_length:,} (1-based)"
        )
    return b1, b2
