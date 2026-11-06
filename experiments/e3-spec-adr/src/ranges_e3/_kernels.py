"""Range math on ``(start, width)`` int64 arrays. No classes (deep module, clear seam).

Every function is pure: arrays in, arrays out. Element indices are 0-based.
"""

from __future__ import annotations

import numpy as np

from ranges_e3._types import BoolArray, IntArray

# --- intra-range (spec 20-intra.md) ---

FIX_START, FIX_END, FIX_CENTER = 0, 1, 2


def shift(start: IntArray, width: IntArray, by: IntArray) -> tuple[IntArray, IntArray]:
    """INTRA.1."""
    return start + by, width


def narrow(
    start: IntArray,
    width: IntArray,
    rel_start: IntArray | None,
    rel_end: IntArray | None,
    rel_width: IntArray | None,
) -> tuple[IntArray, IntArray]:
    """INTRA.2. Relative arguments may be negative; refuses results outside the range."""
    given = (rel_start is not None, rel_end is not None, rel_width is not None)
    if sum(given) == 0 or sum(given) == 3 or given == (False, False, True):
        raise ValueError("narrow needs start and/or end, plus optionally width (not width alone)")
    end = start + width - 1
    a = None if rel_start is None else np.where(rel_start < 0, width + rel_start + 1, rel_start)
    b = None if rel_end is None else np.where(rel_end < 0, width + rel_end + 1, rel_end)
    if a is not None and b is not None:
        new_s, new_e = start + a - 1, start + b - 1
    elif a is not None and rel_width is not None:
        new_s = start + a - 1
        new_e = new_s + rel_width - 1
    elif b is not None and rel_width is not None:
        new_e = start + b - 1
        new_s = new_e - rel_width + 1
    elif a is not None:
        new_s, new_e = start + a - 1, end
    else:
        assert b is not None
        new_s, new_e = start, start + b - 1
    if ((new_s < start) | (new_e > end) | (new_e < new_s - 1)).any():
        raise ValueError("narrow result must lie inside each range")
    return new_s, new_e - new_s + 1


def resize(
    start: IntArray, width: IntArray, new_width: IntArray, fix: IntArray
) -> tuple[IntArray, IntArray]:
    """INTRA.3. ``fix`` holds FIX_* codes per range."""
    if (new_width < 0).any():
        raise ValueError("resize width must be non-negative")
    end = start + width - 1
    new_start = np.where(
        fix == FIX_START,
        start,
        np.where(fix == FIX_END, end - new_width + 1, start + (width - new_width) // 2),
    )
    return new_start, new_width


def flank(
    start: IntArray, width: IntArray, k: IntArray, at_start: BoolArray, both: bool
) -> tuple[IntArray, IntArray]:
    """INTRA.4. ``at_start`` selects the start side per range."""
    end = start + width - 1
    if both:
        a = np.abs(k)
        return np.where(at_start, start - a, end - a + 1), 2 * a
    outward = k >= 0
    new_start = np.where(
        at_start,
        np.where(outward, start - k, start),
        np.where(outward, end + 1, end + k + 1),
    )
    return new_start, np.abs(k)


def promoters(
    start: IntArray, width: IntArray, upstream: IntArray, downstream: IntArray, mirror: BoolArray
) -> tuple[IntArray, IntArray]:
    """INTRA.5. ``mirror`` ranges (``-`` strand) get ``[e-down+1, e+up]``."""
    if (upstream < 0).any() or (downstream < 0).any():
        raise ValueError("promoters upstream/downstream must be non-negative")
    end = start + width - 1
    return np.where(mirror, end - downstream + 1, start - upstream), upstream + downstream


# --- inter-range (spec 30-inter.md, 40-setops.md, 50-coverage.md) ---

_EMPTY = np.empty(0, dtype=np.int64)


def reduce_scan(
    start: IntArray, end: IntArray, min_gap: int, drop_empty: bool
) -> tuple[IntArray, IntArray, IntArray, IntArray]:
    """INTER.2-4. Returns ``(group_start, group_end, order, bounds)``.

    ``order`` lists the kept input indices in scan order ``(start, end, index)``;
    ``bounds`` are the positions in ``order`` where each group begins.
    """
    if min_gap < 0:
        raise ValueError("min_gapwidth must be non-negative")
    n = len(start)
    if n == 0:
        return _EMPTY, _EMPTY, _EMPTY, _EMPTY
    order = np.lexsort((np.arange(n), end, start)).astype(np.int64)
    s, e = start[order], end[order]
    if drop_empty:
        order, s, e, bounds, ge = _scan_dropping_seeds(order, s, e, min_gap)
        return s[bounds] if len(bounds) else _EMPTY, ge, order, bounds
    cummax = np.maximum.accumulate(e)
    new = np.empty(n, dtype=np.bool_)
    new[0] = True
    new[1:] = s[1:] > cummax[:-1] + min_gap
    bounds = np.flatnonzero(new).astype(np.int64)
    last = np.append(bounds[1:] - 1, n - 1)
    return s[bounds], cummax[last], order, bounds


def _scan_dropping_seeds(
    order: IntArray, s: IntArray, e: IntArray, min_gap: int
) -> tuple[IntArray, IntArray, IntArray, IntArray, IntArray]:
    """The sequential scan of INTER.4 (ADR-0007): zero-width ranges may join but not seed."""
    keep: list[int] = []
    bounds: list[int] = []
    ends: list[int] = []
    for i in range(len(s)):
        si, ei = int(s[i]), int(e[i])
        if ends and si <= ends[-1] + min_gap:
            ends[-1] = max(ends[-1], ei)
        elif ei < si:  # zero-width, would seed a group: discard
            continue
        else:
            bounds.append(len(keep))
            ends.append(ei)
        keep.append(i)
    kept = np.array(keep, dtype=np.int64)
    return order[kept], s[kept], e[kept], np.array(bounds, dtype=np.int64), np.array(ends, dtype=np.int64)


def split_groups(order: IntArray, bounds: IntArray) -> list[IntArray]:
    """Per-group index arrays (revmap) from ``reduce_scan`` output."""
    return np.split(order, bounds[1:]) if len(bounds) else []


def span(start: IntArray, end: IntArray) -> tuple[int, int] | None:
    """INTER.1: ``(min start, max end)`` or None when empty."""
    return (int(start.min()), int(end.max())) if len(start) else None


def positive_reduced(start: IntArray, end: IntArray) -> tuple[IntArray, IntArray]:
    """Maximal runs of covered positions (adjacent runs merged); zero-width ranges ignored."""
    pos = end >= start
    gs, ge, _, _ = reduce_scan(start[pos], end[pos], 1, False)
    return gs, ge


def gaps(start: IntArray, end: IntArray, lo: int, hi: int) -> tuple[IntArray, IntArray]:
    """INTER.6: uncovered runs inside ``[lo, hi]``."""
    rs, re = positive_reduced(start, end)
    gs = np.concatenate([[lo], re + 1]).astype(np.int64)
    ge = np.concatenate([rs - 1, [hi]]).astype(np.int64)
    gs, ge = np.maximum(gs, lo), np.minimum(ge, hi)
    ok = gs <= ge
    return gs[ok], ge[ok]


def _piece_cover(
    breaks: IntArray, start: IntArray, end: IntArray
) -> IntArray:
    """Coverage count of each piece ``[breaks[i], breaks[i+1]-1]`` by positive-width ranges."""
    pos = end >= start
    d = np.zeros(len(breaks) + 1, dtype=np.int64)
    np.add.at(d, np.searchsorted(breaks, start[pos]), 1)
    np.add.at(d, np.searchsorted(breaks, end[pos] + 1), -1)
    return np.cumsum(d)[:-1]


def disjoin(
    start: IntArray, end: IntArray, with_revmap: bool
) -> tuple[IntArray, IntArray, list[IntArray]]:
    """INTER.7-8."""
    if len(start) == 0:
        return _EMPTY, _EMPTY, []
    breaks = np.unique(np.concatenate([start, end + 1]))
    cover = _piece_cover(breaks, start, end)[:-1]
    keep = np.flatnonzero(cover > 0)
    ps, pe = breaks[keep], breaks[keep + 1] - 1
    if not with_revmap:
        return ps, pe, []
    pos = np.flatnonzero(end >= start)
    lo = np.searchsorted(breaks, start[pos])
    hi = np.searchsorted(breaks, end[pos] + 1)
    counts = hi - lo
    piece = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)]) if len(pos) else _EMPTY
    owner = np.repeat(pos, counts)
    by_piece = np.argsort(piece, kind="stable")
    piece, owner = piece[by_piece], owner[by_piece]
    revmap = []
    for k in keep:
        a, b = np.searchsorted(piece, k, "left"), np.searchsorted(piece, k, "right")
        revmap.append(owner[a:b].astype(np.int64))
    return ps, pe, revmap


def set_operation(
    xs: IntArray, xe: IntArray, ys: IntArray, ye: IntArray, op: str
) -> tuple[IntArray, IntArray]:
    """SET.1: ``op`` in union/intersect/setdiff over covered positions."""
    rxs, rxe = positive_reduced(xs, xe)
    rys, rye = positive_reduced(ys, ye)
    breaks = np.unique(np.concatenate([rxs, rxe + 1, rys, rye + 1]))
    if len(breaks) < 2:
        return _EMPTY, _EMPTY
    in_x = _piece_cover(breaks, rxs, rxe)[:-1] > 0
    in_y = _piece_cover(breaks, rys, rye)[:-1] > 0
    keep = {"union": in_x | in_y, "intersect": in_x & in_y, "setdiff": in_x & ~in_y}[op]
    prev = np.concatenate([[False], keep[:-1]])
    nxt = np.concatenate([keep[1:], [False]])
    first, last = np.flatnonzero(keep & ~prev), np.flatnonzero(keep & ~nxt)
    return breaks[first], breaks[last + 1] - 1


def punion(
    xs: IntArray, xe: IntArray, ys: IntArray, ye: IntArray, fill_gap: bool
) -> tuple[IntArray, IntArray]:
    """SET.3."""
    gap = np.maximum(xs - ye, ys - xe) - 1
    if not fill_gap and (gap > 0).any():
        raise ValueError("some pair of ranges has a gap; use fill_gap=True")
    return np.minimum(xs, ys), np.maximum(xe, ye)


def pintersect(
    xs: IntArray, xe: IntArray, ys: IntArray, ye: IntArray, resolve_empty: str
) -> tuple[IntArray, IntArray]:
    """SET.4."""
    if resolve_empty not in ("none", "max.start", "start.x"):
        raise ValueError("resolve_empty must be 'none', 'max.start' or 'start.x'")
    s, e = np.maximum(xs, ys), np.minimum(xe, ye)
    bad = e < s - 1
    if bad.any():
        if resolve_empty == "none":
            raise ValueError("some pair of ranges does not intersect; set resolve_empty")
        if resolve_empty == "start.x":
            s = np.where(bad, xs, s)
        e = np.where(bad, s - 1, e)
    return s, e


def psetdiff(
    xs: IntArray, xe: IntArray, ys: IntArray, ye: IntArray
) -> tuple[IntArray, IntArray]:
    """SET.5. Returns ``(start, end)`` per pair; refuses pairs that would split."""
    x_empty, y_empty = xe < xs, ye < ys
    split_by_empty = y_empty & ~x_empty & (xs < ys) & (ys <= xe)
    overlap = ~x_empty & ~y_empty & (np.maximum(xs, ys) <= np.minimum(xe, ye))
    left, right = overlap & (ys > xs), overlap & (ye < xe)
    if (split_by_empty | (left & right)).any():
        raise ValueError("psetdiff would split a range in two")
    s = np.where(right & ~left, ye + 1, xs)
    e = np.where(left & ~right, ys - 1, xe)
    gone = overlap & ~left & ~right
    return s, np.where(gone, xs - 1, e)


def coverage_runs(
    start: IntArray, end: IntArray, shift_by: IntArray, length: int | None
) -> tuple[IntArray, IntArray]:
    """COV.1-3. Run lengths and values of the coverage over positions ``1..L``."""
    s, e = start + shift_by, end + shift_by
    total = length if length is not None else int(max(0, e.max())) if len(e) else 0
    if total == 0:
        return _EMPTY, _EMPTY
    pos = e >= s
    s, e = np.maximum(s[pos], 1), np.minimum(e[pos], total)
    ok = s <= e
    s, e = s[ok], e[ok]
    breaks = np.unique(np.concatenate([[1, total + 1], s, e + 1]))
    values = _piece_cover(breaks, s, e)[:-1]
    lengths = np.diff(breaks)
    change = np.flatnonzero(np.concatenate([[True], values[1:] != values[:-1]]))
    merged = np.add.reduceat(lengths, change) if len(change) else _EMPTY
    return merged.astype(np.int64), values[change].astype(np.int64)
