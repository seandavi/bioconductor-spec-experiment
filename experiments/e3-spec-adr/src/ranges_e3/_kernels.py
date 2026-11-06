"""Range math on ``(start, width)`` int64 arrays. No classes (deep module, clear seam).

Every function is pure: arrays in, arrays out. Element indices are 0-based.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

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
    if sum(given) == 0:
        return start, width
    if sum(given) == 3 or given == (False, False, True):
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
    if min_gap == 0 and n > 1:
        # INTER.3 / ADR-0013: an unabsorbed zero-width range equal to its predecessor vanishes
        twin = new.copy()
        twin[0] = False
        twin[1:] &= (e[1:] < s[1:]) & (e[:-1] < s[:-1]) & (s[1:] == s[:-1])
        if twin.any():
            order, s, e = order[~twin], s[~twin], e[~twin]
            n = len(s)
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


# --- overlaps (spec 60-overlaps.md) ---

OVERLAP_TYPES = ("any", "start", "end", "within", "equal")
_CHUNK = 2_000_000  # max candidate pairs materialised at once


def check_overlap_args(max_gap: int, min_overlap: int, kind: str) -> None:
    """OVL.2-3 argument rules."""
    if kind not in OVERLAP_TYPES:
        raise ValueError(f"type must be one of {OVERLAP_TYPES}")
    if max_gap < -1 or min_overlap < 0:
        raise ValueError("max_gap must be >= -1 and min_overlap >= 0")
    if kind == "any" and max_gap >= 0 and min_overlap > 0:
        raise ValueError("type 'any' accepts max_gap >= 0 or min_overlap > 0, not both")


def find_overlaps(
    qs: IntArray, qe: IntArray, ss: IntArray, se: IntArray, max_gap: int, min_overlap: int, kind: str
) -> tuple[IntArray, IntArray]:
    """OVL.2-7: hit pairs ``(query, subject)`` sorted by query then subject."""
    check_overlap_args(max_gap, min_overlap, kind)
    nq, ns = len(qs), len(ss)
    if nq == 0 or ns == 0:
        return _EMPTY, _EMPTY
    order = np.argsort(ss, kind="stable")
    ss_sorted = ss[order]
    max_width = int((se - ss + 1).max())
    g = max(max_gap, 0)
    if kind == "any":
        lo, hi = qs - max_gap - max_width, qe + max_gap + 1
    elif kind in ("start", "equal"):
        lo, hi = qs - g, qs + g
    elif kind == "end":
        lo, hi = qe - g - max_width + 1, qe + g + 1
    else:  # within
        lo, hi = qe - max_width + 1, qs
    first = np.searchsorted(ss_sorted, lo, "left")
    count = np.maximum(np.searchsorted(ss_sorted, hi, "right") - first, 0)
    ends = np.cumsum(count)
    out_q: list[IntArray] = []
    out_s: list[IntArray] = []
    done = 0
    while done < nq:
        base = ends[done - 1] if done else 0
        stop = int(np.searchsorted(ends, base + _CHUNK, "right"))
        stop = min(max(stop, done + 1), nq)
        q_idx = np.arange(done, stop)
        c = count[done:stop]
        total = int(c.sum())
        if total:
            rep_q = np.repeat(q_idx, c)
            within_q = np.arange(total) - np.repeat(np.cumsum(c) - c, c)
            s_idx = order[np.repeat(first[done:stop], c) + within_q]
            keep = _overlap_condition(
                qs[rep_q], qe[rep_q], ss[s_idx], se[s_idx], max_gap, min_overlap, kind
            )
            out_q.append(rep_q[keep])
            out_s.append(s_idx[keep])
        done = stop
    if not out_q:
        return _EMPTY, _EMPTY
    qi, si = np.concatenate(out_q).astype(np.int64), np.concatenate(out_s).astype(np.int64)
    srt = np.lexsort((si, qi))
    return qi[srt], si[srt]


def _overlap_condition(
    qs: IntArray, qe: IntArray, ss: IntArray, se: IntArray, max_gap: int, min_overlap: int, kind: str
) -> BoolArray:
    g = max(max_gap, 0)
    if kind == "any":
        ok = (np.maximum(qs - se, ss - qe) - 1) <= max_gap
    elif kind == "start":
        ok = np.abs(qs - ss) <= g
    elif kind == "end":
        ok = np.abs(qe - se) <= g
    elif kind == "equal":
        ok = (np.abs(qs - ss) <= g) & (np.abs(qe - se) <= g)
    else:
        ok = (ss <= qs) & (qe <= se)
        if max_gap >= 1:
            ok &= (se - ss) - (qe - qs) <= max_gap
    if min_overlap > 0:
        ok &= np.minimum(qe, se) - np.maximum(qs, ss) + 1 >= min_overlap
    return ok


def masked_ints(values: IntArray, missing: BoolArray) -> np.ma.MaskedArray:  # type: ignore[type-arg]
    """CONV.4: int64 array with NA positions masked."""
    return np.ma.MaskedArray(np.where(missing, 0, values).astype(np.int64), mask=missing)


def select_pairs(qi: IntArray, si: IntArray, nq: int, select: str) -> np.ma.MaskedArray:  # type: ignore[type-arg]
    """OVL.8: per-query first/last/arbitrary subject from pairs sorted by (query, subject)."""
    if select not in ("first", "last", "arbitrary"):
        raise ValueError("select must be 'all', 'first', 'last' or 'arbitrary'")
    out = np.zeros(nq, dtype=np.int64)
    has = np.zeros(nq, dtype=np.bool_)
    if len(qi):
        # numpy fancy assignment keeps the last write per index; order pairs so it is the wanted one
        if select == "last":
            out[qi] = si
        else:
            out[qi[::-1]] = si[::-1]
        has[qi] = True
    return masked_ints(out, ~has)


# --- nearest / precede / follow / distance (spec 70-nearest.md) ---


def distance(qs: IntArray, qe: IntArray, ss: IntArray, se: IntArray) -> IntArray:
    """NEAR.1 (elementwise, arrays already aligned)."""
    return np.maximum(0, np.maximum(qs - se, ss - qe) - 1)


def _expand_groups(
    queries: IntArray, order: IntArray, lo: IntArray, hi: IntArray
) -> tuple[IntArray, IntArray]:
    """Pairs ``(query, order[lo:hi])`` for each listed query."""
    count = hi - lo
    total = int(count.sum())
    if total == 0:
        return _EMPTY, _EMPTY
    within = np.arange(total) - np.repeat(np.cumsum(count) - count, count)
    return np.repeat(queries, count).astype(np.int64), order[np.repeat(lo, count) + within].astype(np.int64)


def precede_pairs(qe: IntArray, ss: IntArray) -> tuple[IntArray, IntArray]:
    """NEAR.3: for each query the subjects sharing the smallest start ``> qe`` (ascending index)."""
    if len(ss) == 0 or len(qe) == 0:
        return _EMPTY, _EMPTY
    order = np.argsort(ss, kind="stable")
    sv = ss[order]
    pos = np.searchsorted(sv, qe, "right")
    ok = np.flatnonzero(pos < len(ss))
    lo = pos[ok]
    return _expand_groups(ok, order, lo, np.searchsorted(sv, sv[lo], "right"))


def follow_pairs(qs: IntArray, se: IntArray) -> tuple[IntArray, IntArray]:
    """NEAR.4: for each query the subjects sharing the largest end ``< qs`` (ascending index)."""
    if len(se) == 0 or len(qs) == 0:
        return _EMPTY, _EMPTY
    order = np.argsort(se, kind="stable")
    ev = se[order]
    pos = np.searchsorted(ev, qs, "left") - 1
    ok = np.flatnonzero(pos >= 0)
    hi = pos[ok] + 1
    return _expand_groups(ok, order, np.searchsorted(ev, ev[pos[ok]], "left"), hi)


def nearest_pairs(
    qs: IntArray, qe: IntArray, ss: IntArray, se: IntArray, ids: tuple[IntArray, IntArray] | None = None
) -> tuple[IntArray, IntArray]:
    """NEAR.2: all subjects at minimal distance, sorted by (query, subject).

    ``ids`` = (query ids, subject ids): pairs with equal ids are excluded (self-search).
    """
    nq, ns = len(qs), len(ss)
    if nq == 0 or ns == 0:
        return _EMPTY, _EMPTY
    q0, s0 = find_overlaps(qs, qe, ss, se, 0, 0, "any")  # every pair at distance 0
    if ids is not None:
        keep = ids[0][q0] != ids[1][s0]
        q0, s0 = q0[keep], s0[keep]
    rest = np.flatnonzero(np.bincount(q0, minlength=nq) == 0)
    by_start, by_end = np.argsort(ss, kind="stable"), np.argsort(se, kind="stable")
    sv, ev = ss[by_start], se[by_end]
    # nearest candidates beyond distance 0: smallest start > qe+1, largest end < qs-1
    right = np.searchsorted(sv, qe[rest] + 1, "right")
    left = np.searchsorted(ev, qs[rest] - 1, "left") - 1
    has_r, has_l = right < ns, left >= 0
    inf = np.iinfo(np.int64).max
    d_r = np.where(has_r, sv[np.minimum(right, ns - 1)] - qe[rest] - 1, inf)
    d_l = np.where(has_l, qs[rest] - ev[np.maximum(left, 0)] - 1, inf)
    best = np.minimum(d_r, d_l)
    take_r = has_r & (d_r == best)
    take_l = has_l & (d_l == best)
    pr = np.flatnonzero(take_r)
    pl_ = np.flatnonzero(take_l)
    rq, rs = _expand_groups(rest[pr], by_start, right[pr], np.searchsorted(sv, sv[right[pr]], "right"))
    lq, ls = _expand_groups(rest[pl_], by_end, np.searchsorted(ev, ev[left[pl_]], "left"), left[pl_] + 1)
    qi = np.concatenate([q0, rq, lq]).astype(np.int64)
    si = np.concatenate([s0, rs, ls]).astype(np.int64)
    srt = np.lexsort((si, qi))
    return qi[srt], si[srt]


# --- ordering (spec 80-order.md) ---


def order_by(keys: list[IntArray], decreasing: bool = False) -> IntArray:
    """ORD.2: indices sorting by ``keys`` (first = most significant); ties keep index order."""
    n = len(keys[0])
    sign = -1 if decreasing else 1
    cols = [sign * k for k in keys]
    return np.lexsort((np.arange(n), *reversed(cols))).astype(np.int64)


def _tie_groups(keys: list[IntArray], order: IntArray) -> IntArray:
    """Group number (0-based, along ``order``) of each sorted element; equal keys share one."""
    n = len(order)
    new = np.zeros(n, dtype=np.bool_)
    if n:
        new[0] = True
        for k in keys:
            ks = k[order]
            new[1:] |= ks[1:] != ks[:-1]
    return np.cumsum(new) - 1


def rank_by(keys: list[IntArray], method: str) -> IntArray | npt.NDArray[np.float64]:
    """ORD.4: 1-based ranks."""
    if method not in ("average", "min", "max", "first"):
        raise ValueError("ties_method must be 'average', 'min', 'max' or 'first'")
    n = len(keys[0])
    order = order_by(keys)
    group = _tie_groups(keys, order)
    position = np.arange(n, dtype=np.int64)
    first_pos = np.zeros(n, dtype=np.int64)
    last_pos = np.zeros(n, dtype=np.int64)
    if n:
        starts = np.flatnonzero(np.concatenate([[True], group[1:] != group[:-1]]))
        ends = np.append(starts[1:], n) - 1
        first_pos, last_pos = starts[group], ends[group]
    if method == "first":
        sorted_rank: npt.NDArray[np.generic] = position + 1
    elif method == "min":
        sorted_rank = first_pos + 1
    elif method == "max":
        sorted_rank = last_pos + 1
    else:
        sorted_rank = (first_pos + last_pos) / 2 + 1
    out = np.empty(n, dtype=sorted_rank.dtype)
    out[order] = sorted_rank
    return out  # type: ignore[return-value]


def duplicated_by(keys: list[IntArray], from_last: bool) -> BoolArray:
    """ORD.5."""
    n = len(keys[0])
    order = order_by(keys)
    group = _tie_groups(keys, order)
    dup_sorted = np.zeros(n, dtype=np.bool_)
    if n:
        same_as_prev = np.concatenate([[False], group[1:] == group[:-1]])
        same_as_next = np.concatenate([group[1:] == group[:-1], [False]])
        dup_sorted = same_as_next if from_last else same_as_prev
    out = np.empty(n, dtype=np.bool_)
    out[order] = dup_sorted
    return out
