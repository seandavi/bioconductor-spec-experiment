"""Ordering kernels on a tuple of integer key arrays (primary key first)."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def order_k(keys: Sequence[np.ndarray], decreasing: bool = False) -> np.ndarray:
    """Stable ordering by ``keys`` (first key most significant).

    ``decreasing=True`` reverses the key order but keeps tied elements in their
    original relative order, like R's ``order(decreasing=TRUE)``.
    """
    n = len(keys[0])
    if n == 0:
        return np.array([], dtype=np.int64)
    ks = [-k if decreasing else k for k in keys]
    return np.lexsort(tuple(reversed(ks))).astype(np.int64)


def _tie_groups(keys: Sequence[np.ndarray], order: np.ndarray) -> np.ndarray:
    """For the sorted sequence, the group number (0-based) of each position."""
    n = len(order)
    if n == 0:
        return np.array([], dtype=np.int64)
    new = np.zeros(n, dtype=bool)
    new[0] = True
    for k in keys:
        ks = k[order]
        new[1:] |= ks[1:] != ks[:-1]
    return np.cumsum(new) - 1


def rank_k(keys: Sequence[np.ndarray], ties_method: str = "first") -> np.ndarray:
    """1-based ranks of the elements (ascending). Returns float for ``"average"``."""
    n = len(keys[0])
    if ties_method not in ("average", "first", "last", "min", "max"):
        raise ValueError("'ties_method' must be one of average, first, last, min, max")
    if n == 0:
        return np.array([], dtype=float if ties_method == "average" else np.int64)
    if ties_method == "last":
        o = order_k([k[::-1] for k in keys])
        o = (n - 1 - o).astype(np.int64)  # positions in the original array
        pos = np.empty(n, dtype=np.int64)
        pos[o] = np.arange(1, n + 1)
        return pos
    o = order_k(keys)
    grp = _tie_groups(keys, o)
    sorted_rank = np.arange(1, n + 1)
    if ties_method == "first":
        out_sorted = sorted_rank
    else:
        counts = np.bincount(grp)
        ends = np.cumsum(counts)  # last 1-based rank of each group
        firsts = ends - counts + 1
        if ties_method == "min":
            out_sorted = firsts[grp]
        elif ties_method == "max":
            out_sorted = ends[grp]
        else:
            out_sorted = ((firsts + ends) / 2)[grp]
    out = np.empty(n, dtype=out_sorted.dtype)
    out[o] = out_sorted
    return out


def duplicated_k(keys: Sequence[np.ndarray], from_last: bool = False) -> np.ndarray:
    """True for every element equal (on all keys) to an earlier one
    (or, with ``from_last``, to a later one)."""
    n = len(keys[0])
    if n == 0:
        return np.array([], dtype=bool)
    if from_last:
        return duplicated_k([k[::-1] for k in keys])[::-1].copy()
    o = order_k(keys)  # stable: within ties original order is kept
    grp = _tie_groups(keys, o)
    first_in_group = np.concatenate(([True], grp[1:] != grp[:-1]))
    dup = np.empty(n, dtype=bool)
    dup[o] = ~first_in_group
    return dup
