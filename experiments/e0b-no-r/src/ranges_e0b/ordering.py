"""Ordering kernels over lists of integer key columns (most significant first)."""

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

from ._sew import BoolArray, IntArray, empty_int


def order_keys(keys: Sequence[IntArray], decreasing: bool = False) -> IntArray:
    """Stable ordering permutation (0-based); ties keep input order, also if decreasing."""
    if len(keys[0]) == 0:
        return empty_int()
    cols = [-k if decreasing else k for k in keys]
    return np.lexsort(cols[::-1]).astype(np.int64)


def _tie_groups(keys: Sequence[IntArray], order: IntArray) -> IntArray:
    """Tie-group id (0..g-1) of each element in sorted order."""
    n = len(order)
    new = np.zeros(n, dtype=bool)
    if n:
        new[0] = True
        for k in keys:
            ks = k[order]
            new[1:] |= ks[1:] != ks[:-1]
    return np.cumsum(new) - 1


def rank_keys(keys: Sequence[IntArray], method: str = "average") -> npt.NDArray[np.float64] | IntArray:
    """1-based ranks. "average" returns float64; others int64."""
    if method not in ("average", "first", "min", "max"):
        raise ValueError(f"unknown ties_method {method!r}")
    n = len(keys[0])
    order = order_keys(keys)
    gid = _tie_groups(keys, order)
    pos = np.arange(1, n + 1, dtype=np.int64)
    if n:
        first = np.zeros(gid[-1] + 1, dtype=np.int64)
        last = np.zeros_like(first)
        # sorted positions are increasing, so assignment order picks min / max
        first[gid[::-1]] = pos[::-1]
        last[gid] = pos
    out: npt.NDArray[np.float64] | IntArray
    if method == "average":
        out = np.empty(n, dtype=np.float64)
        out[order] = (first[gid] + last[gid]) / 2 if n else []
    else:
        out = np.empty(n, dtype=np.int64)
        src = {"first": pos, "min": first[gid] if n else pos, "max": last[gid] if n else pos}[method]
        out[order] = src
    return out


def duplicated_keys(keys: Sequence[IntArray]) -> BoolArray:
    """True for every element equal to an earlier one (first occurrence is False)."""
    n = len(keys[0])
    dup = np.zeros(n, dtype=bool)
    if n == 0:
        return dup
    order = order_keys(keys)
    gid = _tie_groups(keys, order)
    same_as_prev = gid[1:] == gid[:-1]
    dup[order[1:][same_as_prev]] = True
    return dup
