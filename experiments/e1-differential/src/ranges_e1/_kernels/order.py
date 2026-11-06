"""Ordering, ranking and duplicate detection over lexicographic integer keys."""

from __future__ import annotations

import numpy as np

TIES = ("average", "first", "last", "random", "max", "min")


def order(keys: list[np.ndarray], decreasing: bool = False) -> np.ndarray:
    """Stable ordering by keys (most significant first); ties keep index order either way."""
    n = len(keys[0])
    if n == 0:
        return np.empty(0, dtype=np.int64)
    ks = [-k if decreasing else k for k in keys]
    return np.lexsort(ks[::-1]).astype(np.int64)


def _tie_groups(keys: list[np.ndarray], ordr: np.ndarray) -> np.ndarray:
    """Group id (0..g-1) of each element in sorted order."""
    if len(ordr) == 0:
        return np.empty(0, dtype=np.int64)
    change = np.zeros(len(ordr), dtype=bool)
    change[0] = True
    for k in keys:
        s = k[ordr]
        change[1:] |= s[1:] != s[:-1]
    return np.cumsum(change) - 1


def rank(keys: list[np.ndarray], ties_method: str) -> np.ndarray:
    """R's rank(): 1-based ranks (float for "average", int otherwise)."""
    n = len(keys[0])
    ordr = order(keys)
    grp = _tie_groups(keys, ordr)
    pos = np.arange(1, n + 1)
    first_pos = np.zeros(n, dtype=np.int64)
    last_pos = np.zeros(n, dtype=np.int64)
    if n:
        starts = np.flatnonzero(np.concatenate([[True], grp[1:] != grp[:-1]]))
        ends = np.append(starts[1:], n) - 1
        first_pos = (starts + 1)[grp]
        last_pos = (ends + 1)[grp]
    if ties_method == "min":
        out_sorted = first_pos
    elif ties_method == "max":
        out_sorted = last_pos
    elif ties_method == "first":
        out_sorted = pos
    elif ties_method == "last":
        out_sorted = first_pos + last_pos - pos
    elif ties_method == "average":
        out_sorted = (first_pos + last_pos) / 2
    else:
        raise ValueError(f"'arg' should be one of {', '.join(TIES)}")
    out = np.empty(n, dtype=out_sorted.dtype)
    out[ordr] = out_sorted
    return out


def duplicated(keys: list[np.ndarray], from_last: bool = False) -> np.ndarray:
    n = len(keys[0])
    if n == 0:
        return np.zeros(0, dtype=bool)
    ordr = order(keys)
    grp = _tie_groups(keys, ordr)  # ties keep index order inside a group
    dup_sorted = np.zeros(n, dtype=bool)
    if from_last:
        dup_sorted[:-1] = grp[:-1] == grp[1:]
    else:
        dup_sorted[1:] = grp[1:] == grp[:-1]
    out = np.zeros(n, dtype=bool)
    out[ordr] = dup_sorted
    return out
