"""Dependency-free columnar cache for canonical annotation frames.

Stores a DataFrame as a compressed ``.npz``: numeric/bool columns as raw
arrays, categorical columns as integer codes plus a UTF-8 category table.
Loads in O(seconds) for multi-million-row genomes without requiring
pyarrow/parquet. The format is versioned; a mismatch triggers re-parsing
from the retained source file rather than a hard failure.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

CACHE_FORMAT_VERSION = 1


class FrameCacheError(RuntimeError):
    pass


def save_frame(df: pd.DataFrame, path: str | os.PathLike) -> None:
    arrays: dict[str, np.ndarray] = {
        "__format_version__": np.array([CACHE_FORMAT_VERSION], dtype="int64"),
        "__n_rows__": np.array([len(df)], dtype="int64"),
        "__columns__": np.array(list(df.columns), dtype="U64"),
    }
    for col in df.columns:
        s = df[col]
        if isinstance(s.dtype, pd.CategoricalDtype):
            arrays[f"cat_codes::{col}"] = s.cat.codes.to_numpy(dtype="int32")
            arrays[f"cat_values::{col}"] = np.asarray(s.cat.categories, dtype="U")
        elif s.dtype == object:
            arrays[f"str::{col}"] = s.to_numpy(dtype="U")
        else:
            arrays[f"num::{col}"] = s.to_numpy()
    tmp = os.fspath(path) + ".tmp"
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, **arrays)
    os.replace(tmp, path)


def load_frame(path: str | os.PathLike) -> pd.DataFrame:
    try:
        npz = np.load(path, allow_pickle=False)
    except (ValueError, OSError, EOFError) as e:
        raise FrameCacheError(f"{os.fspath(path)}: unreadable frame cache ({e})") from e
    with npz as z:
        try:
            version = int(z["__format_version__"][0])
        except KeyError:
            raise FrameCacheError(f"{os.fspath(path)}: not a TEscope frame cache") from None
        if version != CACHE_FORMAT_VERSION:
            raise FrameCacheError(
                f"{os.fspath(path)}: cache format v{version} != expected v{CACHE_FORMAT_VERSION}"
            )
        columns = [str(c) for c in z["__columns__"]]
        data = {}
        for col in columns:
            if f"num::{col}" in z:
                data[col] = z[f"num::{col}"]
            elif f"cat_codes::{col}" in z:
                codes = z[f"cat_codes::{col}"]
                cats = [str(v) for v in z[f"cat_values::{col}"]]
                data[col] = pd.Categorical.from_codes(codes, categories=cats)
            elif f"str::{col}" in z:
                data[col] = [str(v) for v in z[f"str::{col}"]]
            else:
                raise FrameCacheError(f"{os.fspath(path)}: missing column payload for {col!r}")
        return pd.DataFrame(data, columns=columns)
