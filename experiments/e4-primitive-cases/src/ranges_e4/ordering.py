"""Ordering kernels over a key matrix ``keys`` of shape (k, n), most significant key first.

IRanges key = (start, width); GRanges key = (seqlevel code, strand code, start, width).
All sorts are stable, so ties keep input order (also when decreasing, as in R).
"""

from __future__ import annotations

import numpy as np

TIES = ("average", "first", "last", "min", "max")


def order(keys: np.ndarray, decreasing: bool = False) -> np.ndarray:
    if keys.shape[1] == 0:
        return np.zeros(0, np.int64)
    k = -keys if decreasing else keys
    return np.lexsort(k[::-1]).astype(np.int64)


def _groups(keys: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(order, group start flags in sorted order) for stable ascending sort."""
    o = order(keys)
    sk = keys[:, o]
    new = np.r_[True, (sk[:, 1:] != sk[:, :-1]).any(axis=0)] if len(o) else np.zeros(0, bool)
    return o, new


def duplicated(keys: np.ndarray, from_last: bool = False) -> np.ndarray:
    n = keys.shape[1]
    o, new = _groups(keys)
    out = np.zeros(n, dtype=bool)
    if from_last:
        # last of each group (in index order) is the original
        is_last = np.r_[new[1:], True] if n else new
        out[o[~is_last]] = True
    else:
        out[o[~new]] = True
    return out


def rank(keys: np.ndarray, ties_method: str) -> np.ndarray:
    if ties_method not in TIES:
        raise ValueError(f"'ties_method' must be one of {TIES}")
    n = keys.shape[1]
    o, new = _groups(keys)
    gid = np.cumsum(new) - 1
    first = np.flatnonzero(new)
    sizes = np.diff(np.r_[first, n])
    lo = first[gid] + 1  # min rank of the group (1-based), in sorted order
    hi = lo + sizes[gid] - 1
    pos = np.arange(n)
    out = np.empty(n, dtype=np.float64 if ties_method == "average" else np.int64)
    match ties_method:
        case "min":
            r = lo
        case "max":
            r = hi
        case "average":
            r = (lo + hi) / 2
        case "first":
            r = pos + 1
        case _:  # last: later occurrences get lower ranks
            r = lo + hi - (pos + 1)
    out[o] = r
    return out
