"""Subscript resolution shared by IRanges and GRanges (CORE.9)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from ranges_e3._types import IntArray


def resolve(index: Any, n: int, names: Sequence[str | None] | None) -> IntArray:
    """Translate any supported subscript into a validated int64 position array."""
    if isinstance(index, slice):
        return np.arange(n, dtype=np.int64)[index]
    if isinstance(index, (int, np.integer)):
        return _positions(np.array([index], dtype=np.int64), n)
    if isinstance(index, str):
        return _by_name([index], names)
    if isinstance(index, np.ndarray):
        arr = index
    else:
        seq = list(index)
        if seq and all(isinstance(v, str) for v in seq):
            return _by_name(seq, names)
        arr = np.asarray(seq)
    if arr.size == 0:
        return np.empty(0, dtype=np.int64)
    if arr.dtype == np.bool_:
        if arr.shape != (n,):
            raise IndexError(f"boolean mask of length {arr.size}, expected {n}")
        return np.flatnonzero(arr).astype(np.int64)
    if arr.dtype.kind in "US":
        return _by_name([str(v) for v in arr], names)
    if arr.dtype.kind not in "iu":
        raise IndexError("subscript must be integers, booleans or names")
    return _positions(arr.astype(np.int64), n)


def _positions(pos: IntArray, n: int) -> IntArray:
    if pos.size and (pos.min() < -n or pos.max() >= n):
        raise IndexError(f"position out of range for length {n}")
    return np.where(pos < 0, pos + n, pos)


def _by_name(wanted: Sequence[str], names: Sequence[str | None] | None) -> IntArray:
    if names is None:
        raise IndexError("cannot subscript by name: object has no names")
    first: dict[str, int] = {}
    for i, name in enumerate(names):
        if name is not None and name not in first:
            first[name] = i
    missing = [w for w in wanted if w not in first]
    if missing:
        raise IndexError(f"unknown names: {missing[:3]}")
    return np.array([first[w] for w in wanted], dtype=np.int64)
