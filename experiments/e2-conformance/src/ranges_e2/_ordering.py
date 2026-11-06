"""Ordering kernels over lists of int64 key arrays (most significant key first)."""

from __future__ import annotations

import numpy as np

from ._util import IntArray

TIES_METHODS = ("average", "first", "last", "min", "max")


def order_by(keys: list[IntArray], decreasing: bool = False) -> IntArray:
    """Stable ordering; with decreasing=True keys are reversed but ties keep input order."""
    if not keys or len(keys[0]) == 0:
        return np.zeros(0, dtype=np.int64)
    ks = [-k for k in keys] if decreasing else keys
    return np.lexsort(tuple(ks[::-1])).astype(np.int64)


def rank_by(keys: list[IntArray], ties_method: str = "average") -> np.ndarray:
    """R's rank over composite keys. average -> float64; others -> int64 (1-based)."""
    if ties_method not in TIES_METHODS:
        raise ValueError(f"ties_method must be one of {TIES_METHODS}")
    n = len(keys[0]) if keys else 0
    order = order_by(keys)
    if n == 0:
        return np.zeros(0, dtype=np.float64 if ties_method == "average" else np.int64)
    stacked = np.stack([k[order] for k in keys])
    new_group = np.ones(n, dtype=bool)
    new_group[1:] = np.any(stacked[:, 1:] != stacked[:, :-1], axis=0)
    group_id = np.cumsum(new_group) - 1
    first = np.nonzero(new_group)[0]
    sizes = np.diff(np.concatenate((first, [n])))
    lo = first[group_id] + 1  # min rank per element (1-based)
    hi = lo + sizes[group_id] - 1
    pos = np.arange(n) - first[group_id]  # position within tie group
    if ties_method == "min":
        sorted_rank = lo
    elif ties_method == "max":
        sorted_rank = hi
    elif ties_method == "first":
        sorted_rank = lo + pos
    elif ties_method == "last":
        sorted_rank = hi - pos
    else:
        sorted_rank = (lo + hi) / 2.0
    out = np.empty(n, dtype=sorted_rank.dtype)
    out[order] = sorted_rank
    return out


def duplicated_by(keys: list[IntArray], from_last: bool = False) -> np.ndarray:
    n = len(keys[0]) if keys else 0
    if n == 0:
        return np.zeros(0, dtype=bool)
    stacked = np.stack(keys, axis=1)
    if from_last:
        stacked = stacked[::-1]
    _, first_idx = np.unique(stacked, axis=0, return_index=True)
    dup = np.ones(n, dtype=bool)
    dup[first_idx] = False
    return dup[::-1] if from_last else dup
