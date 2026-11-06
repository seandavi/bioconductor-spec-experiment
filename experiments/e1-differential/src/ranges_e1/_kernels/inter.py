"""Inter-range ops (range / reduce / gaps / disjoin) on flat arrays with a group key.

Every function takes `group` (int64, any order) alongside `start`, `width`;
results come back sorted by (group, start). IRanges passes a single group.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_SHIFT = 1 << 33  # separates groups when taking running maxima across them


@dataclass(frozen=True)
class Reduced:
    group: np.ndarray
    start: np.ndarray
    width: np.ndarray
    members: list[np.ndarray]  # input positions per output range, in sweep order


def range_(group: np.ndarray, start: np.ndarray, width: np.ndarray) -> Reduced:
    """Min start / max end per group (empty ranges included)."""
    if len(start) == 0:
        return Reduced(group, start, width, [])
    end = start + width - 1
    order = np.lexsort((start, group))
    g = group[order]
    first = np.concatenate([[True], g[1:] != g[:-1]])
    gids = g[first]
    bounds = np.flatnonzero(first)
    mn = np.minimum.reduceat(start[order], bounds)
    mx = np.maximum.reduceat(end[order], bounds)
    members = [order[b:e] for b, e in zip(bounds, np.append(bounds[1:], len(g)))]
    return Reduced(gids, mn, mx - mn + 1, members)


def reduce(
    group: np.ndarray,
    start: np.ndarray,
    width: np.ndarray,
    min_gapwidth: int = 1,
    drop_empty: bool = False,
) -> Reduced:
    """Merge ranges whose gap is `< min_gapwidth`; stable sweep by (group, start, width)."""
    pos = np.arange(len(start), dtype=np.int64)
    if drop_empty:
        keep = width > 0
        group, start, width, pos = group[keep], start[keep], width[keep], pos[keep]
    if len(start) == 0:
        e = np.empty(0, dtype=np.int64)
        return Reduced(e, e, e, [])
    order = np.lexsort((width, start, group))
    g, s, w, p = group[order], start[order], width[order], pos[order]
    end = s + w - 1
    run_end = np.maximum.accumulate(g * _SHIFT + end + (_SHIFT >> 1)) - g * _SHIFT - (_SHIFT >> 1)
    new = np.ones(len(s), dtype=bool)
    new[1:] = (g[1:] != g[:-1]) | (s[1:] - run_end[:-1] - 1 >= min_gapwidth)
    firsts = np.flatnonzero(new)
    lasts = np.append(firsts[1:], len(s)) - 1
    out_start = s[firsts]
    out_end = run_end[lasts]
    out_g, out_w = g[firsts], out_end - out_start + 1
    members = [p[a : b + 1] for a, b in zip(firsts, lasts)]
    # R quirk: with min_gapwidth=0 identical empty ranges form separate clusters that
    # R then drops as repeats of the previous output range (their members vanish too).
    dup = np.zeros(len(out_g), dtype=bool)
    dup[1:] = (out_g[1:] == out_g[:-1]) & (out_start[1:] == out_start[:-1]) & (out_w[1:] == out_w[:-1])
    if dup.any():
        keep = ~dup
        members = [m for m, k in zip(members, keep) if k]
        out_g, out_start, out_w = out_g[keep], out_start[keep], out_w[keep]
    return Reduced(out_g, out_start, out_w, members)


def gaps(
    group: np.ndarray,
    start: np.ndarray,
    width: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Uncovered stretches inside [lo[g], hi[g]] for every group id 0..len(lo)-1.

    Empty ranges cover nothing and are dropped first (as R does).
    """
    red = reduce(group, start, width, 1, drop_empty=True)
    out_g: list[np.ndarray] = []
    out_s: list[np.ndarray] = []
    out_w: list[np.ndarray] = []
    for gid in range(len(lo)):
        sel = red.group == gid
        s = red.start[sel]
        e = s + red.width[sel] - 1
        a = np.concatenate([[lo[gid]], e + 1])
        b = np.concatenate([s - 1, [hi[gid]]])
        a = np.maximum(a, lo[gid])
        b = np.minimum(b, hi[gid])
        ok = a <= b
        out_g.append(np.full(int(ok.sum()), gid, dtype=np.int64))
        out_s.append(a[ok])
        out_w.append(b[ok] - a[ok] + 1)
    cat = lambda xs: np.concatenate(xs) if xs else np.empty(0, dtype=np.int64)  # noqa: E731
    return cat(out_g), cat(out_s), cat(out_w)


def disjoin(
    group: np.ndarray, start: np.ndarray, width: np.ndarray, with_members: bool = False
) -> Reduced:
    """Segments between all start / end+1 breakpoints that some non-empty range covers.

    Empty ranges add breakpoints but are never reported.
    """
    e = np.empty(0, dtype=np.int64)
    if len(start) == 0:
        return Reduced(e, e, e, [])
    end1 = start + width  # exclusive end
    out_g: list[int] = []
    out_s: list[int] = []
    out_w: list[int] = []
    members: list[np.ndarray] = []
    for gid in np.unique(group):
        idx = np.flatnonzero(group == gid)
        s, e1, w = start[idx], end1[idx], width[idx]
        bps = np.unique(np.concatenate([s, e1]))
        if len(bps) < 2:
            continue
        ne = w > 0
        delta = np.zeros(len(bps), dtype=np.int64)
        np.add.at(delta, np.searchsorted(bps, s[ne]), 1)
        np.add.at(delta, np.searchsorted(bps, e1[ne]), -1)
        cover = np.cumsum(delta)[:-1]
        for k in np.flatnonzero(cover > 0):
            out_g.append(int(gid))
            out_s.append(int(bps[k]))
            out_w.append(int(bps[k + 1] - bps[k]))
            if with_members:
                members.append(idx[ne & (s <= bps[k]) & (e1 >= bps[k + 1])])
    return Reduced(
        np.array(out_g, dtype=np.int64), np.array(out_s, dtype=np.int64), np.array(out_w, dtype=np.int64), members
    )
