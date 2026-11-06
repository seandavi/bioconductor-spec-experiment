"""Inter-range kernels on (start, width) int64 arrays: reduce, gaps, disjoin, set operations.

Ranges are closed intervals [start, start + width - 1]; zero-width ranges have end = start - 1.
All functions return NEW arrays; indices in revmaps are 0-based positions into the inputs.
"""

from __future__ import annotations

import numpy as np

from ._util import IntArray, cummax

Revmap = list[IntArray]


def _empty() -> IntArray:
    return np.zeros(0, dtype=np.int64)


def reduce_ranges(
    start: IntArray,
    width: IntArray,
    min_gapwidth: int = 1,
    drop_empty: bool = False,
    with_revmap: bool = False,
) -> tuple[IntArray, IntArray, Revmap | None]:
    """Merge ranges whose gap is < min_gapwidth. Output sorted by start.

    A range starts a new group when `start - (max end+1 so far) >= min_gapwidth`; identical
    zero-width ranges always collapse (R quirk, pinned by IR.REDUCE.EMPTY).
    revmap members are listed in (start, width, index) order.
    """
    if min_gapwidth < 0:
        raise ValueError("min.gapwidth must be non-negative")
    idx = np.arange(len(start), dtype=np.int64)
    if drop_empty:
        keep = width > 0
        idx = idx[keep]
    if len(idx) == 0:
        return _empty(), _empty(), ([] if with_revmap else None)
    order = idx[np.lexsort((width[idx], start[idx]))]
    s = start[order]
    w = width[order]
    excl = s + w
    prev_max = np.concatenate(([s[0]], cummax(excl)[:-1]))
    gap = s - prev_max
    new_group = gap >= min_gapwidth
    same_empty = np.zeros(len(s), dtype=bool)
    same_empty[1:] = (w[1:] == 0) & (w[:-1] == 0) & (s[1:] == s[:-1]) & (gap[1:] == 0)
    new_group &= ~same_empty
    new_group[0] = True
    first = np.nonzero(new_group)[0]
    out_start = s[first]
    out_excl = np.maximum.reduceat(excl, first)
    revmap = None
    if with_revmap:
        revmap = [chunk for chunk in np.split(order, first[1:])]
    return out_start, out_excl - out_start, revmap


def gaps_ranges(
    start: IntArray, width: IntArray, lo: int | None, hi: int | None
) -> tuple[IntArray, IntArray]:
    """Uncovered positions within [lo, hi]. Defaults: first start / last end of the non-empty ranges."""
    rs, rw, _ = reduce_ranges(start, width, 1, drop_empty=True)
    if lo is None:
        if len(rs) == 0:
            return _empty(), _empty()
        lo = int(rs[0])
    if hi is None:
        if len(rs) == 0:
            return _empty(), _empty()
        hi = int(rs[-1] + rw[-1] - 1)
    gs = np.concatenate(([lo], rs + rw)).astype(np.int64)
    ge = np.concatenate((rs - 1, [hi])).astype(np.int64)
    gs = np.maximum(gs, lo)
    ge = np.minimum(ge, hi)
    keep = ge >= gs
    return gs[keep], ge[keep] - gs[keep] + 1


def disjoin_ranges(
    start: IntArray, width: IntArray, with_revmap: bool = False
) -> tuple[IntArray, IntArray, Revmap | None]:
    """Split into non-overlapping pieces (adjacent pieces stay separate); zero-width ranges only split."""
    keep = np.nonzero(width > 0)[0]
    if len(keep) == 0:
        return _empty(), _empty(), ([] if with_revmap else None)
    s = start[keep]
    e1 = s + width[keep]  # exclusive end
    # zero-width ranges add no coverage but still split the pieces they sit inside (R behaviour)
    points = np.unique(np.concatenate((s, e1, start[width == 0])))
    first = np.searchsorted(points, s)
    last = np.searchsorted(points, e1)  # interval k is [points[k], points[k+1])
    delta = np.zeros(len(points) + 1, dtype=np.int64)
    np.add.at(delta, first, 1)
    np.add.at(delta, last, -1)
    covered = np.cumsum(delta)[: len(points) - 1] > 0
    ks = np.nonzero(covered)[0]
    out_start = points[ks]
    out_width = points[ks + 1] - points[ks]
    revmap = None
    if with_revmap:
        members: list[list[int]] = [[] for _ in range(len(points) - 1)]
        for r in range(len(keep)):
            for k in range(first[r], last[r]):
                members[k].append(int(keep[r]))
        revmap = [np.array(members[k], dtype=np.int64) for k in ks]
    return out_start, out_width, revmap


def _coverage_pair(
    xs: IntArray, xw: IntArray, ys: IntArray, yw: IntArray
) -> tuple[IntArray, IntArray, IntArray]:
    """Elementary intervals from the union of breakpoints and per-set coverage flags."""
    points = np.unique(np.concatenate((xs, xs + xw, ys, ys + yw)))
    if len(points) < 2:
        return points, np.zeros(0, dtype=bool), np.zeros(0, dtype=bool)

    def covered(s: IntArray, w: IntArray) -> IntArray:
        delta = np.zeros(len(points) + 1, dtype=np.int64)
        if len(s):
            np.add.at(delta, np.searchsorted(points, s), 1)
            np.add.at(delta, np.searchsorted(points, s + w), -1)
        return np.cumsum(delta)[: len(points) - 1] > 0

    return points, covered(xs, xw), covered(ys, yw)


def _intervals_from_mask(points: IntArray, mask: IntArray) -> tuple[IntArray, IntArray]:
    ks = np.nonzero(mask)[0]
    if len(ks) == 0:
        return _empty(), _empty()
    return reduce_ranges(points[ks], points[ks + 1] - points[ks], 1)[:2]


def intersect_ranges(xs: IntArray, xw: IntArray, ys: IntArray, yw: IntArray) -> tuple[IntArray, IntArray]:
    xs, xw, _ = reduce_ranges(xs, xw, 1, drop_empty=True)
    ys, yw, _ = reduce_ranges(ys, yw, 1, drop_empty=True)
    if len(xs) == 0 or len(ys) == 0:
        return _empty(), _empty()
    points, cx, cy = _coverage_pair(xs, xw, ys, yw)
    return _intervals_from_mask(points, cx & cy)


def setdiff_ranges(xs: IntArray, xw: IntArray, ys: IntArray, yw: IntArray) -> tuple[IntArray, IntArray]:
    xs, xw, _ = reduce_ranges(xs, xw, 1, drop_empty=True)
    ys, yw, _ = reduce_ranges(ys, yw, 1, drop_empty=True)
    if len(xs) == 0:
        return _empty(), _empty()
    if len(ys) == 0:
        return xs, xw
    points, cx, cy = _coverage_pair(xs, xw, ys, yw)
    return _intervals_from_mask(points, cx & ~cy)


def union_ranges(xs: IntArray, xw: IntArray, ys: IntArray, yw: IntArray) -> tuple[IntArray, IntArray]:
    s = np.concatenate((xs, ys))
    w = np.concatenate((xw, yw))
    return reduce_ranges(s, w, 1, drop_empty=True)[:2]


def range_of(start: IntArray, width: IntArray) -> tuple[IntArray, IntArray]:
    """Single range from min start to max end (zero-width ranges count with end = start - 1)."""
    if len(start) == 0:
        return _empty(), _empty()
    lo = start.min()
    hi = (start + width - 1).max()
    return np.array([lo], dtype=np.int64), np.array([hi - lo + 1], dtype=np.int64)
