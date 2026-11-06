"""Shared helpers: int32 coercion, R-style recycling, subscript normalisation."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt

IntArray = npt.NDArray[np.int64]

INT32_MAX = 2**31 - 1


def as_int64(values: Any, what: str = "value") -> IntArray:
    """Coerce numbers to int64. Floats truncate toward zero (R's as.integer); NaN/inf/NA error."""
    if isinstance(values, np.ndarray) and values.dtype.kind in "iu":
        return values.astype(np.int64)
    if isinstance(values, (int, float, np.integer, np.floating)):
        values = [values]
    arr = np.asarray(values)
    if arr.size == 0:
        return np.zeros(0, dtype=np.int64)
    if arr.dtype == object:
        if any(v is None for v in arr.tolist()):
            raise ValueError(f"{what} cannot contain NAs")
        arr = np.asarray(arr.tolist())
    if arr.dtype.kind == "b":
        return arr.astype(np.int64)
    if arr.dtype.kind == "f":
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{what} cannot contain NA/NaN/Inf")
        arr = np.trunc(arr)
    elif arr.dtype.kind not in "iu":
        raise TypeError(f"{what} must be numeric, got {arr.dtype}")
    return arr.astype(np.int64)


def check_int32(a: IntArray, what: str = "value") -> IntArray:
    if a.size and (a.max() > INT32_MAX or a.min() < -INT32_MAX):
        raise ValueError(f"{what} overflows 32-bit integer range")
    return a


def rep_len(a: np.ndarray, n: int, what: str = "argument") -> np.ndarray:
    """R's rep_len: cyclic recycling to length n. Empty only recycles to empty."""
    if len(a) == n:
        return a
    if len(a) == 0:
        raise ValueError(f"{what} has length 0 but {n} values are needed")
    return np.resize(a, n)


def as_bool_array(values: Any, n: int, what: str = "argument") -> npt.NDArray[np.bool_]:
    arr = np.asarray(values)
    if arr.dtype.kind != "b":
        raise TypeError(f"{what} must be logical")
    return rep_len(arr, n, what)


def positions_from_subscript(i: Any, n: int, names: Sequence[str | None] | None) -> IntArray:
    """Normalise a Python-API subscript to 0-based positions (numpy semantics).

    int / slice / bool array (length n) / int array (negative counts from end) / list of str.
    """
    if isinstance(i, slice):
        return np.arange(n, dtype=np.int64)[i]
    if isinstance(i, (int, np.integer)) and not isinstance(i, (bool, np.bool_)):
        k = int(i)
        if k < -n or k >= n:
            raise IndexError(f"index {k} out of range for length {n}")
        return np.array([k % n if n else k], dtype=np.int64)
    if isinstance(i, (str, bytes)):
        i = [i]
    if isinstance(i, (list, tuple)) and len(i) > 0 and all(isinstance(e, str) for e in i):
        if names is None:
            raise KeyError("cannot subset by name: object has no names")
        lookup: dict[str, int] = {}
        for pos, nm in enumerate(names):
            if nm is not None and nm not in lookup:
                lookup[nm] = pos
        try:
            return np.array([lookup[e] for e in i], dtype=np.int64)
        except KeyError as e:
            raise KeyError(f"unknown name {e.args[0]!r}") from None
    arr = np.asarray(i)
    if arr.size == 0:
        return np.zeros(0, dtype=np.int64)
    if arr.dtype.kind == "b":
        if len(arr) != n:
            raise IndexError(f"logical subscript has length {len(arr)}, expected {n}")
        return np.nonzero(arr)[0].astype(np.int64)
    if arr.dtype.kind in "iu":
        arr = arr.astype(np.int64)
        if arr.min() < -n or arr.max() >= n:
            raise IndexError("subscript contains out-of-bounds indices")
        return np.where(arr < 0, arr + n, arr)
    raise TypeError(f"unsupported subscript type {arr.dtype}")


def cummax(a: IntArray) -> IntArray:
    return np.maximum.accumulate(a) if len(a) else a
