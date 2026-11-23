"""Display-only signal transforms: smoothing and per-copy quantiles.

These functions NEVER touch canonical results. Smoothing is a
visualization transform applied at render/derived-export time; the
unsmoothed profile remains the scientific output. Quantiles summarize the
existing per-copy projection matrix (variation ACROSS TE COPIES — not the
statistical uncertainty of the mean, which is the SEM).

All methods are NaN-aware: bins without data stay NaN and never borrow
values from neighbors' absence (normalized convolution for MA/Gaussian;
Savitzky-Golay fits only within all-finite windows, pure numpy).
"""

from __future__ import annotations

import numpy as np

SMOOTHING_METHODS = ("none", "moving_average", "gaussian", "savgol")


class TransformError(ValueError):
    pass


def validate_smoothing(method: str, window: int, sigma: float, polyorder: int) -> None:
    if method not in SMOOTHING_METHODS:
        raise TransformError(f"smoothing method must be one of {SMOOTHING_METHODS}")
    if method == "moving_average" and (window < 3 or window % 2 == 0):
        raise TransformError("moving-average window must be an odd integer >= 3")
    if method == "gaussian" and sigma <= 0:
        raise TransformError("gaussian sigma must be > 0")
    if method == "savgol":
        if window < 5 or window % 2 == 0:
            raise TransformError("Savitzky-Golay window must be an odd integer >= 5")
        if not (1 <= polyorder < window):
            raise TransformError("Savitzky-Golay polyorder must satisfy 1 <= order < window")


def smooth(values: np.ndarray, method: str, window: int = 5,
           sigma: float = 2.0, polyorder: int = 2) -> np.ndarray:
    """Return the smoothed copy of ``values`` (float64); input untouched."""
    validate_smoothing(method, window, sigma, polyorder)
    v = np.asarray(values, dtype="float64").copy()
    if method == "none" or len(v) == 0:
        return v
    if method == "moving_average":
        kernel = np.ones(window)
    elif method == "gaussian":
        half = max(int(np.ceil(3 * sigma)), 1)
        xs = np.arange(-half, half + 1)
        kernel = np.exp(-0.5 * (xs / sigma) ** 2)
    else:  # savgol
        return _savgol_nan(v, window, polyorder)
    return _nan_convolve(v, kernel)


def _same_convolve(a: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """np.convolve 'same' that always returns len(a) (numpy returns
    max(len(a), len(kernel)) when the kernel is longer)."""
    res = np.convolve(a, kernel, mode="same")
    if len(res) > len(a):
        off = (len(res) - len(a)) // 2
        res = res[off:off + len(a)]
    return res


def _nan_convolve(v: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    finite = np.isfinite(v)
    filled = np.where(finite, v, 0.0)
    num = _same_convolve(filled, kernel)
    den = _same_convolve(finite.astype("float64"), kernel)
    out = np.full_like(v, np.nan)
    ok = den > 1e-12
    out[ok] = num[ok] / den[ok]
    out[~finite] = np.nan  # never invent values where the profile has none
    return out


def _savgol_nan(v: np.ndarray, window: int, polyorder: int) -> np.ndarray:
    """Savitzky-Golay via precomputed least-squares coefficients (numpy only).

    The center-point coefficient vector c solves the classic SG design:
    smoothed[i] = c . v[i-h : i+h+1]. Windows containing any NaN fall back
    to the NaN-aware moving average for that point (documented behavior);
    NaN bins stay NaN.
    """
    half = window // 2
    x = np.arange(-half, half + 1, dtype="float64")
    A = np.vander(x, polyorder + 1, increasing=True)  # (window, order+1)
    # projection row for the window center: e0^T (A^T A)^-1 A^T
    coeffs = np.linalg.pinv(A)[0]  # evaluates fitted polynomial at x=0
    n = len(v)
    finite = np.isfinite(v)
    out = np.full_like(v, np.nan)
    ma = _nan_convolve(v, np.ones(window))
    for i in range(n):
        if not finite[i]:
            continue
        lo, hi = i - half, i + half + 1
        if lo < 0 or hi > n or not finite[lo:hi].all():
            out[i] = ma[i]
        else:
            out[i] = float(coeffs @ v[lo:hi])
    return out


def copy_quantiles(matrix, qs: tuple[float, ...]) -> dict[float, np.ndarray]:
    """Per-bin quantiles of per-copy values from a CopyBinMatrix.

    Summarizes variation across TE copies (NaN where <2 observations).
    """
    vals = np.asarray(matrix.values, dtype="float64")
    n_obs = np.isfinite(vals).sum(axis=0)
    out: dict[float, np.ndarray] = {}
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for q in qs:
            arr = np.nanpercentile(vals, q, axis=0)
            arr[n_obs < 2] = np.nan
            out[q] = arr
    return out
