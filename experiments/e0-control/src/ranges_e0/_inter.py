"""Inter-range kernels (range / reduce / gaps / disjoin) on ``(start, width)`` arrays."""

from __future__ import annotations

import numpy as np

I64 = np.int64


def range_k(start: np.ndarray, width: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One range ``[min(start), max(end)]``; nothing for empty input."""
    if len(start) == 0:
        return start.copy(), width.copy()
    s = start.min()
    e = (start + width - 1).max()
    return np.array([s], dtype=I64), np.array([e - s + 1], dtype=I64)


def reduce_k(
    start: np.ndarray,
    width: np.ndarray,
    min_gapwidth: int = 1,
    drop_empty: bool = False,
    with_revmap: bool = False,
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray] | None]:
    """Merge ranges separated by a gap narrower than ``min_gapwidth``.

    Sweep over the ranges sorted by (start, end): a range joins the current
    group when ``start <= running_max_end + min_gapwidth``. Empty ranges take
    part in the sweep unless ``drop_empty`` (so an isolated empty range
    survives and one inside/adjacent to a real range is absorbed, like R).
    Returns ``(start, width, revmap)``; ``revmap`` (0-based input indices per
    output range, in sweep order) only when requested.
    """
    idx = np.arange(len(start), dtype=I64)
    if drop_empty:
        idx = idx[width > 0]
    if len(idx) == 0:
        return (
            np.array([], dtype=I64),
            np.array([], dtype=I64),
            [] if with_revmap else None,
        )
    s_all = start[idx]
    e_all = s_all + width[idx] - 1
    o = np.lexsort((idx, e_all, s_all))  # start, then end, then index
    s = s_all[o]
    e = e_all[o]
    ce = np.maximum.accumulate(e)
    new_group = np.empty(len(s), dtype=bool)
    new_group[0] = True
    new_group[1:] = s[1:] > ce[:-1] + min_gapwidth
    first = np.flatnonzero(new_group)
    last = np.concatenate((first[1:] - 1, [len(s) - 1]))
    out_start = s[first]
    out_width = ce[last] - out_start + 1
    revmap = None
    if with_revmap:
        members = idx[o]
        revmap = np.split(members, first[1:])
    return out_start.astype(I64), out_width.astype(I64), revmap


def gaps_k(
    start: np.ndarray,
    width: np.ndarray,
    gap_start: int | None = None,
    gap_end: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Complement of the union of the ranges within ``[gap_start, gap_end]``.

    Missing bounds default to the first start / last end of the non-empty
    ranges; with no non-empty range and a missing bound there are no gaps (R)."""
    rs, rw, _ = reduce_k(start, width, 1, drop_empty=True)
    re_ = rs + rw - 1
    if len(rs) == 0 and (gap_start is None or gap_end is None):
        return np.array([], dtype=I64), np.array([], dtype=I64)
    gs = int(rs[0]) if gap_start is None else int(gap_start)
    ge = int(re_[-1]) if gap_end is None else int(gap_end)
    g_start = np.concatenate(([gs], re_ + 1)).astype(I64)
    g_end = np.concatenate((rs - 1, [ge])).astype(I64)
    g_start = np.maximum(g_start, gs)
    g_end = np.minimum(g_end, ge)
    keep = g_end >= g_start
    return g_start[keep], (g_end - g_start + 1)[keep]


def disjoin_k(
    start: np.ndarray, width: np.ndarray, with_revmap: bool = False
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray] | None]:
    """Disjoint pieces of the covered positions, split at every range boundary.

    Breakpoints are all starts and all ``end + 1``; empty ranges contribute a
    breakpoint (at their start) but no coverage, which is what R does."""
    n = len(start)
    if n == 0:
        return start.copy(), width.copy(), [] if with_revmap else None
    end1 = start + width  # end + 1
    bp = np.unique(np.concatenate((start, end1)))
    seg_cover = np.zeros(len(bp), dtype=I64)
    ne = width > 0
    ne_idx = np.flatnonzero(ne)
    lo = np.searchsorted(bp, start[ne_idx])
    hi = np.searchsorted(bp, end1[ne_idx])
    np.add.at(seg_cover, lo, 1)
    np.add.at(seg_cover, hi, -1)
    cover = np.cumsum(seg_cover)[:-1]  # segment i = [bp[i], bp[i+1] - 1]
    keep = np.flatnonzero(cover > 0)
    out_start = bp[keep]
    out_width = bp[keep + 1] - bp[keep]
    revmap = None
    if with_revmap:
        counts = hi - lo  # segments spanned by each non-empty range
        seg_of_pair = np.repeat(lo, counts) + (
            np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
        )
        range_of_pair = np.repeat(ne_idx, counts)
        o = np.lexsort((range_of_pair, seg_of_pair))
        seg_sorted = seg_of_pair[o]
        range_sorted = range_of_pair[o]
        cuts = np.flatnonzero(np.diff(seg_sorted)) + 1
        groups = np.split(range_sorted, cuts)
        seg_ids = seg_sorted[np.concatenate(([0], cuts))] if len(seg_sorted) else np.array([], dtype=I64)
        by_seg = dict(zip(seg_ids.tolist(), groups))
        revmap = [by_seg[i] for i in keep.tolist()]
    return out_start.astype(I64), out_width.astype(I64), revmap
