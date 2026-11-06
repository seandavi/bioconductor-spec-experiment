"""Inter-range kernels (range, reduce, gaps, disjoin) over start/width arrays.

Reduce is one sorted sweep. Ranges are visited in ``order()`` (start, then
width, then index); a range joins the current run when
``start - run_end - 1 < min_gapwidth``. Zero-width ranges follow the same
rule, which reproduces R: an empty range touching a non-empty one merges at
the default ``min_gapwidth=1`` but not at 0.
"""

from __future__ import annotations

import numpy as np

EMPTY = np.zeros(0, dtype=np.int64)


def sort_order(start: np.ndarray, width: np.ndarray) -> np.ndarray:
    """Positions ordering by (start, width), stable on ties."""
    return np.lexsort((width, start))


def range_of(start: np.ndarray, width: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if not len(start):
        return EMPTY, EMPTY
    lo = start.min()
    hi = (start + width - 1).max()
    return np.array([lo]), np.array([hi - lo + 1])


def reduce_runs(
    start: np.ndarray, width: np.ndarray, min_gapwidth: int, drop_empty: bool
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """Return (run_start, run_width, members); members index the input arrays."""
    if min_gapwidth < 0:
        raise ValueError("'min_gapwidth' must be non-negative")
    keep = np.flatnonzero(width > 0) if drop_empty else np.arange(len(start))
    if not len(keep):
        return EMPTY, EMPTY, []
    keep = keep[sort_order(start[keep], width[keep])]
    s = start[keep]
    e = s + width[keep] - 1
    cum_end = np.maximum.accumulate(e)
    new_run = np.concatenate([[True], s[1:] - cum_end[:-1] - 1 >= min_gapwidth])
    first = np.flatnonzero(new_run)
    last = np.concatenate([first[1:], [len(s)]]) - 1
    run_start = s[first]
    run_width = cum_end[last] - run_start + 1
    members = np.split(keep, first[1:])
    return run_start, run_width, members


def gaps(
    start: np.ndarray, width: np.ndarray, lo: int | None, hi: int | None
) -> tuple[np.ndarray, np.ndarray]:
    """Complement of the non-empty coverage within [lo, hi].

    A missing bound defaults to the coverage extent; with no coverage and a
    missing bound the result is empty (as in R).
    """
    rs, rw, _ = reduce_runs(start, width, 1, True)
    re = rs + rw - 1
    if lo is None and len(rs):
        lo = int(rs[0])
    if hi is None and len(rs):
        hi = int(re[-1])
    if lo is None or hi is None:
        return EMPTY, EMPTY
    g_start = np.concatenate([[lo], re + 1])
    g_end = np.concatenate([rs - 1, [hi]])
    g_start = np.maximum(g_start, lo)
    g_end = np.minimum(g_end, hi)
    ok = g_end >= g_start
    return g_start[ok], g_end[ok] - g_start[ok] + 1


def disjoin(
    start: np.ndarray, width: np.ndarray
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """Maximal segments with a constant, non-empty covering set; ranges of width 0 are ignored."""
    idx = np.flatnonzero(width > 0)
    if not len(idx):
        return EMPTY, EMPTY, []
    s = start[idx]
    e = s + width[idx]  # exclusive end
    bp = np.unique(np.concatenate([s, e]))
    lo = np.searchsorted(bp, s)
    hi = np.searchsorted(bp, e)
    delta = np.zeros(len(bp) + 1, dtype=np.int64)
    np.add.at(delta, lo, 1)
    np.add.at(delta, hi, -1)
    covered = np.cumsum(delta)[:-1] > 0
    seg = np.flatnonzero(covered)
    seg = seg[seg < len(bp) - 1]
    members: list[list[int]] = [[] for _ in range(len(bp) - 1)]
    for pos, (a, b) in enumerate(zip(lo.tolist(), hi.tolist(), strict=True)):
        for k in range(a, b):
            members[k].append(int(idx[pos]))
    out_start = bp[seg]
    out_width = bp[seg + 1] - bp[seg]
    return out_start, out_width, [np.array(members[k], dtype=np.int64) for k in seg.tolist()]
