"""Argument coercion helpers shared by IRanges and GRanges."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
import polars as pl


def as_int_array(x: Any, name: str = "argument") -> np.ndarray:
    """Coerce a scalar / sequence / array to a 1-d int64 array.

    Integral floats are accepted (R accepts doubles); non-integral floats and
    NaN raise ``ValueError``.
    """
    arr = np.asarray(x)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise ValueError(f"'{name}' must be one-dimensional")
    if arr.dtype == np.bool_:
        raise TypeError(f"'{name}' must be integer, not bool")
    if arr.dtype.kind == "f":
        if not np.all(np.isfinite(arr)) or not np.all(arr == np.floor(arr)):
            raise ValueError(f"'{name}' must contain whole numbers (no NaN/fractions)")
        return arr.astype(np.int64)
    if arr.dtype.kind in "iu":
        return arr.astype(np.int64, copy=False)
    if arr.size == 0:
        return arr.astype(np.int64)
    raise TypeError(f"'{name}' must be integer-valued, got dtype {arr.dtype}")


def recycle(arr: np.ndarray, n: int, name: str) -> np.ndarray:
    """Recycle a length-1 array to length ``n``; otherwise lengths must match."""
    if len(arr) == n:
        return arr
    if len(arr) == 1:
        return np.repeat(arr, n)
    raise ValueError(f"'{name}' has length {len(arr)} but {n} was expected")


def as_names(names: Iterable[str | None] | None, n: int) -> np.ndarray | None:
    """Coerce names to an object array of str (``None`` entries become ``""``)."""
    if names is None:
        return None
    arr = np.empty(n, dtype=object)
    vals = list(names)
    if len(vals) != n:
        raise ValueError(f"'names' has length {len(vals)} but {n} was expected")
    for i, v in enumerate(vals):
        arr[i] = "" if v is None else str(v)
    return arr


def as_mcols(mcols: pl.DataFrame | dict[str, Sequence[Any]] | None, n: int) -> pl.DataFrame:
    """Coerce metadata columns to a polars DataFrame of height ``n``."""
    if mcols is None:
        return _empty_mcols(n)
    df = mcols if isinstance(mcols, pl.DataFrame) else pl.DataFrame(mcols)
    if df.width == 0:
        return _empty_mcols(n)
    if df.height != n:
        raise ValueError(f"'mcols' has {df.height} rows but object has length {n}")
    return df


def _empty_mcols(n: int) -> pl.DataFrame:
    """Zero-column mcols (polars frames with no columns have no height)."""
    return pl.DataFrame()


def take_rows(mcols: pl.DataFrame, idx: np.ndarray) -> pl.DataFrame:
    """Select rows of an mcols frame (no-op for zero columns)."""
    if mcols.width == 0:
        return mcols
    return mcols[idx]


def concat_mcols(frames: Sequence[pl.DataFrame], lengths: Sequence[int]) -> pl.DataFrame:
    """Row-bind mcols aligned by column name; missing columns become null."""
    total = int(sum(lengths))
    cols: list[str] = []
    for f in frames:
        for c in f.columns:
            if c not in cols:
                cols.append(c)
    if not cols:
        return _empty_mcols(total)
    parts = []
    for f, n in zip(frames, lengths):
        if f.width == 0:
            f = pl.DataFrame({c: pl.Series(c, [None] * n) for c in cols})
        parts.append(f)
    return pl.concat(parts, how="diagonal_relaxed").select(cols)


def resolve_index(
    index: Any, n: int, names: np.ndarray | None
) -> np.ndarray:
    """Turn an int / slice / int array / bool array / str / list[str] index into positions."""
    if isinstance(index, slice):
        return np.arange(n, dtype=np.int64)[index]
    if isinstance(index, (int, np.integer)) and not isinstance(index, bool):
        i = int(index)
        if not -n <= i < n:
            raise IndexError(f"index {i} out of range for length {n}")
        return np.array([i % n], dtype=np.int64) if n else np.array([], dtype=np.int64)
    if isinstance(index, str):
        index = [index]
    arr = np.asarray(index)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.size == 0:
        return np.array([], dtype=np.int64)
    if arr.dtype == np.bool_:
        if len(arr) != n:
            raise IndexError(f"boolean index has length {len(arr)} but object has length {n}")
        return np.flatnonzero(arr).astype(np.int64)
    if arr.dtype.kind in "OUS":
        if names is None:
            raise KeyError("cannot subset by name when names are None")
        lookup: dict[str, int] = {}
        for i, nm in enumerate(names):
            if nm != "" and nm not in lookup:
                lookup[nm] = i
        out = np.empty(len(arr), dtype=np.int64)
        for j, key in enumerate(arr):
            if key not in lookup:
                raise KeyError(f"subscript contains invalid name: {key!r}")
            out[j] = lookup[key]
        return out
    if arr.dtype.kind in "iu":
        idx = arr.astype(np.int64)
        if np.any((idx < -n) | (idx >= n)):
            raise IndexError("subscript contains out-of-bounds indices")
        return np.where(idx < 0, idx + n, idx)
    raise TypeError(f"unsupported index type {arr.dtype}")
