"""precede / follow / nearest / distance kernels (formulae follow R's IRanges source)."""

from __future__ import annotations

import numpy as np

from ._overlaps import overlap_pairs, select_hits

I64 = np.int64
NO_MATCH = -1
Pairs = tuple[np.ndarray, np.ndarray]


def _expand(qidx: np.ndarray, lo: np.ndarray, hi: np.ndarray, order: np.ndarray) -> Pairs:
    """Pairs ``(q, order[k])`` for every ``k`` in ``[lo, hi)`` of each query in ``qidx``."""
    counts = hi - lo
    total = int(counts.sum())
    if total == 0:
        return np.array([], dtype=I64), np.array([], dtype=I64)
    q = np.repeat(qidx, counts)
    first = np.cumsum(counts) - counts
    k = np.arange(total, dtype=I64) - np.repeat(first, counts) + np.repeat(lo, counts)
    return q.astype(I64), order[k].astype(I64)


def precede_k(xs: np.ndarray, xw: np.ndarray, ss: np.ndarray, select: str = "first") -> np.ndarray | Pairs:
    """Subject range(s) with the smallest start strictly greater than ``end(x)``.

    ``first``: lowest subject index among ties on that start; ``all``: every tie."""
    n = len(xs)
    xe = xs + xw - 1
    order = np.argsort(ss, kind="stable")
    s = ss[order]
    i = np.searchsorted(s, xe, side="right")
    valid = i < len(s)
    if select == "first":
        out = np.full(n, NO_MATCH, dtype=I64)
        out[valid] = order[i[valid]]
        return out
    q = np.flatnonzero(valid)
    iv = i[valid]
    hi = np.searchsorted(s, s[iv], side="right") if len(iv) else iv
    return _expand(q, iv, hi, order)


def follow_k(xs: np.ndarray, ss: np.ndarray, sw: np.ndarray, select: str = "last") -> np.ndarray | Pairs:
    """Subject range(s) with the largest end strictly smaller than ``start(x)``.

    ``last``: highest subject index among ties on that end; ``all``: every tie."""
    n = len(xs)
    se = ss + sw - 1
    order = np.argsort(se, kind="stable")
    e = se[order]
    j = np.searchsorted(e, xs - 1, side="right") - 1
    valid = j >= 0
    if select == "last":
        out = np.full(n, NO_MATCH, dtype=I64)
        out[valid] = order[j[valid]]
        return out
    q = np.flatnonzero(valid)
    jv = j[valid]
    lo = np.searchsorted(e, e[jv], side="left") if len(jv) else jv
    return _expand(q, lo, jv + 1, order)


def nearest_k(
    xs: np.ndarray,
    xw: np.ndarray,
    ss: np.ndarray,
    sw: np.ndarray,
    select: str = "arbitrary",
    drop_self: bool = False,
) -> np.ndarray | Pairs:
    """R's ``nearest``: any overlapping-or-adjacent subject wins (distance 0);
    otherwise the closer of ``precede`` / ``follow``. Ties between the two sides
    take the ``follow`` side for ``arbitrary`` and both sides for ``all``.
    ``arbitrary`` picks the first overlapping hit."""
    n = len(xs)
    if len(ss) == 0:
        if select == "all":
            return np.array([], dtype=I64), np.array([], dtype=I64)
        return np.full(n, NO_MATCH, dtype=I64)
    qi, si = overlap_pairs(xs, xw, ss, sw, "any", 0, 0)
    if drop_self:
        keep = qi != si
        qi, si = qi[keep], si[keep]
    has_hit = np.zeros(n, dtype=bool)
    has_hit[qi] = True
    rest = np.flatnonzero(~has_hit)
    xs_r, xw_r = xs[rest], xw[rest]
    xe_r = xs_r + xw_r - 1
    se = ss + sw - 1
    before = precede_k(xs_r, xw_r, ss, "first")
    after = follow_k(xs_r, ss, sw, "last")
    assert isinstance(before, np.ndarray) and isinstance(after, np.ndarray)
    b_ok, a_ok = before >= 0, after >= 0
    big = np.iinfo(I64).max // 4
    leftdist = np.where(b_ok, ss[np.where(b_ok, before, 0)] - xe_r, big)
    rightdist = np.where(a_ok, xs_r - se[np.where(a_ok, after, 0)], big)
    both = b_ok & a_ok
    left = np.where(both, leftdist < rightdist, ~a_ok)
    if select == "arbitrary":
        out = select_hits(qi, si, n, "first", NO_MATCH)
        out[rest] = np.where(left, before, after)
        return out
    right = ~left
    left = np.where(both & (leftdist == rightdist), True, left)
    ball = precede_k(xs_r, xw_r, ss, "all")
    aall = follow_k(xs_r, ss, sw, "all")
    assert isinstance(ball, tuple) and isinstance(aall, tuple)
    bq, bs = ball
    aq, as_ = aall
    bkeep = left[bq] if len(bq) else np.array([], dtype=bool)
    akeep = right[aq] if len(aq) else np.array([], dtype=bool)
    # bq / aq index into `rest`
    q_all = np.concatenate((qi, rest[bq[bkeep]], rest[aq[akeep]]))
    s_all = np.concatenate((si, bs[bkeep], as_[akeep]))
    o = np.lexsort((s_all, q_all))
    return q_all[o], s_all[o]


def distance_k(xs: np.ndarray, xw: np.ndarray, ys: np.ndarray, yw: np.ndarray) -> np.ndarray:
    """R >= IRanges 2.47.3 ``distance``: 0 when the ranges overlap, else gap width + 1
    (adjacent ranges give 1). Empty ranges never overlap (R tests ``max(start) <= min(end)``)."""
    if len(xs) == 0 or len(ys) == 0:
        return np.array([], dtype=I64)
    n = max(len(xs), len(ys))
    for a, b in ((len(xs), "x"), (len(ys), "y")):
        if a not in (1, n):
            raise ValueError(f"'{b}' must have length 1 or the length of the other argument")
    max_start = np.maximum(xs, ys)
    min_end = np.minimum(xs + xw - 1, ys + yw - 1)
    return np.maximum(max_start - min_end, 0).astype(I64)
