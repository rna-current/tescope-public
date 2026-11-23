"""bigWig continuous tracks with a dual-backend reader.

``pyBigWig`` is the battle-tested C reader but publishes manylinux wheels
only; ``pybigtools`` (Rust bigtools bindings) ships wheels for
Linux/macOS/Windows including arm64. TEscope prefers pyBigWig when
importable (Linux installs) and falls back to pybigtools, behind one
small interface. Golden tests run against every importable backend and
assert identical numbers, so results do not depend on the platform.
"""

from __future__ import annotations

import os

import numpy as np

from .base import ContinuousTrack, TrackKind, TrackMeta

_BACKENDS: list[str] = []
try:  # pragma: no cover - import probing
    import pyBigWig  # type: ignore

    _BACKENDS.append("pybigwig")
except ImportError:  # pragma: no cover
    pyBigWig = None
try:  # pragma: no cover
    import pybigtools  # type: ignore

    _BACKENDS.append("pybigtools")
except ImportError:  # pragma: no cover
    pybigtools = None


def available_backends() -> list[str]:
    return list(_BACKENDS)


class BigWigReadError(RuntimeError):
    pass


class _PyBigWigReader:
    def __init__(self, path: str):
        self._bw = pyBigWig.open(path)
        if self._bw is None or not self._bw.isBigWig():
            raise BigWigReadError(f"{path}: not a bigWig file")

    def chroms(self) -> dict[str, int]:
        return dict(self._bw.chroms())

    def values(self, chrom: str, start: int, end: int) -> np.ndarray:
        clen = self._bw.chroms(chrom)
        out = np.full(end - start, np.nan, dtype="float64")
        if not clen:
            return out
        a, b = max(0, start), min(int(clen), end)
        if b > a:
            vals = self._bw.values(chrom, a, b, numpy=True)
            out[a - start : b - start] = vals
        return out

    def close(self):
        self._bw.close()


class _PyBigToolsReader:
    def __init__(self, path: str):
        # Open via a file object rather than the path string: pybigtools'
        # Rust URL detection misreads Windows drive-letter paths ("C:\...")
        # as URL schemes. File-object reads are numerically identical
        # (asserted by the backend-parity tests).
        self._fh = None
        try:
            self._fh = open(path, "rb")
            self._bw = pybigtools.open(self._fh)
        except Exception as e:  # noqa: BLE001 - backend raises plain Exception
            if self._fh is not None:
                self._fh.close()
                self._fh = None
            # fallback for pybigtools versions without file-object support
            try:
                self._bw = pybigtools.open(path)
            except Exception:  # noqa: BLE001
                raise BigWigReadError(f"{path}: cannot open as bigWig ({e})") from e
        self._chroms = {str(k): int(v) for k, v in self._bw.chroms().items()}

    def chroms(self) -> dict[str, int]:
        return dict(self._chroms)

    def values(self, chrom: str, start: int, end: int) -> np.ndarray:
        out = np.full(end - start, np.nan, dtype="float64")
        clen = self._chroms.get(chrom)
        if clen is None:
            return out
        a, b = max(0, start), min(clen, end)
        if b > a:
            try:
                vals = self._bw.values(chrom, a, b, fillna=np.nan)
            except TypeError:  # pybigtools < 0.3 spelled the kwarg 'missing'
                vals = self._bw.values(chrom, a, b, missing=np.nan)
            out[a - start : b - start] = np.asarray(vals, dtype="float64")
        return out

    def close(self):
        self._bw.close()
        if self._fh is not None:
            self._fh.close()


def open_bigwig_reader(path: str, backend: str | None = None):
    """Open a low-level reader; backend None = first available."""
    if backend is None:
        if not _BACKENDS:
            raise BigWigReadError(
                "no bigWig backend available (install pyBigWig or pybigtools)"
            )
        backend = _BACKENDS[0]
    if backend == "pybigwig":
        if pyBigWig is None:
            raise BigWigReadError("pyBigWig not installed")
        return _PyBigWigReader(path)
    if backend == "pybigtools":
        if pybigtools is None:
            raise BigWigReadError("pybigtools not installed")
        return _PyBigToolsReader(path)
    raise ValueError(f"unknown bigWig backend {backend!r}")


class BigWigTrack(ContinuousTrack):
    def __init__(self, meta: TrackMeta, backend: str | None = None):
        super().__init__(meta)
        if not os.path.exists(meta.path):
            raise FileNotFoundError(meta.path)
        self._reader = open_bigwig_reader(meta.path, backend)
        self.backend = type(self._reader).__name__

    def chroms(self) -> dict[str, int]:
        return self._reader.chroms()

    def values(self, chrom: str, start: int, end: int) -> np.ndarray:
        return self._reader.values(chrom, start, end)

    def close(self) -> None:
        self._reader.close()


def open_bigwig(path: str, name: str | None = None, **meta_kwargs) -> BigWigTrack:
    meta = TrackMeta(
        name=name or os.path.basename(path),
        path=os.fspath(path),
        kind=TrackKind.CONTINUOUS,
        format="bigwig",
        **meta_kwargs,
    )
    return BigWigTrack(meta)
