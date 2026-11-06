"""Small array helpers shared by the range classes."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

INT32_MAX = 2**31 - 1
INT32_MIN = -(2**31) + 1  # R integer NA is -2^31


def as_int(values: Any, what: str) -> np.ndarray:
    """Integer vector (int64) that fits R's int32; None/NaN is rejected."""
    if isinstance(values, np.ma.MaskedArray):
        if np.ma.getmaskarray(values).any():
            raise ValueError(f"{what}: NA not allowed")
        values = values.data
    if isinstance(values, (int, np.integer)):
        values = [values]
    arr = np.asarray(values)
    if arr.size == 0:
        return np.zeros(0, dtype=np.int64)
    if arr.dtype == object or arr.dtype.kind == "f":
        if any(v is None or (isinstance(v, float) and v != v) for v in arr.tolist()):
            raise ValueError(f"{what}: NA not allowed")
        if arr.dtype.kind == "f" and not np.all(arr == np.floor(arr)):
            raise ValueError(f"{what}: non-integer value")
        arr = arr.astype(np.int64)
    if arr.dtype.kind not in "iub":
        raise TypeError(f"{what}: expected integers, got {arr.dtype}")
    arr = arr.astype(np.int64)
    if arr.size and (arr.max() > INT32_MAX or arr.min() < INT32_MIN):
        raise ValueError(f"{what}: value outside R integer range")
    return arr


def recycle(arrs: Sequence[np.ndarray | None]) -> list[np.ndarray | None]:
    """Recycle length-1 arrays to the common length; other mismatches raise."""
    n = max((len(a) for a in arrs if a is not None), default=0)
    out: list[np.ndarray | None] = []
    for a in arrs:
        if a is None or len(a) == n:
            out.append(a)
        elif len(a) == 1:
            out.append(np.repeat(a, n))
        elif n and len(a) == 0:
            raise ValueError("zero-length argument cannot be recycled")
        else:
            raise ValueError(f"argument lengths differ ({len(a)} vs {n})")
    return out
