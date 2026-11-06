"""Inter-range kernels on grouped (start, end, group) arrays.

`g` is an int group code (always 0 for IRanges; seqlevel*3 + strand for
GRanges). Groups never interact. Results are sorted by (group, start).
"""

import numpy as np

from ._sew import IntArray, empty_int

Triple = tuple[IntArray, IntArray, IntArray]


def _empty3() -> Triple:
    return empty_int(), empty_int(), empty_int()


def reduce_kernel(
    s: IntArray,
    e: IntArray,
    g: IntArray,
    min_gapwidth: int = 1,
    drop_empty: bool = False,
    with_revmap: bool = False,
) -> tuple[IntArray, IntArray, IntArray, list[IntArray] | None]:
    """Merge ranges whose gap (next start - running max end - 1) < min_gapwidth."""
    if min_gapwidth < 0:
        raise ValueError("min_gapwidth must be non-negative")
    n = len(s)
    if n == 0:
        return (*_empty3(), [] if with_revmap else None)
    order = np.lexsort((e, s, g))
    ss, ee, gg = s[order], e[order], g[order]
    same_group = gg[1:] == gg[:-1]
    grank = np.concatenate(([0], np.cumsum(~same_group)))
    emin = ee.min()
    span = ee.max() - emin + 1
    # offsetting by group rank makes one global running max behave per group
    cm = np.maximum.accumulate((ee - emin) + grank * span)
    prev_end = cm[:-1] - grank[1:] * span + emin
    new_run = np.ones(n, dtype=bool)
    new_run[1:] = ~same_group | (ss[1:] - prev_end - 1 >= min_gapwidth)
    first = np.flatnonzero(new_run)
    last = np.append(first[1:] - 1, n - 1)
    out_s, out_g = ss[first], gg[first]
    out_e = cm[last] - grank[last] * span + emin
    revmap: list[IntArray] | None = None
    if with_revmap:
        run_id = np.cumsum(new_run) - 1
        by_run = order[np.lexsort((order, run_id))]  # ascending index within run
        revmap = np.split(by_run, first[1:])
    if drop_empty:
        keep = out_e >= out_s
        out_s, out_e, out_g = out_s[keep], out_e[keep], out_g[keep]
        if revmap is not None:
            revmap = [r for r, k in zip(revmap, keep) if k]
    return out_s, out_e, out_g, revmap


def range_kernel(s: IntArray, e: IntArray, g: IntArray) -> Triple:
    """Per group [min start, max end] (empty ranges included)."""
    if len(s) == 0:
        return _empty3()
    o = np.argsort(g, kind="stable")
    gs = g[o]
    b = np.flatnonzero(np.concatenate(([True], gs[1:] != gs[:-1])))
    return np.minimum.reduceat(s[o], b), np.maximum.reduceat(e[o], b), gs[b]


def _unique_points(g: IntArray, v: IntArray) -> tuple[IntArray, IntArray, IntArray, IntArray]:
    """Sort (g, v) points; return order, unique g, unique v, and the unique id of each sorted point."""
    order = np.lexsort((v, g))
    gs, vs = g[order], v[order]
    new = np.ones(len(gs), dtype=bool)
    new[1:] = (gs[1:] != gs[:-1]) | (vs[1:] != vs[:-1])
    uid = np.cumsum(new) - 1
    return order, gs[new], vs[new], uid


def disjoin_kernel(s: IntArray, e: IntArray, g: IntArray) -> Triple:
    """Disjoint pieces, split at every start and end+1; only covered pieces kept.

    Empty ranges add breakpoints (like R) but are never output.
    """
    n = len(s)
    if n == 0:
        return _empty3()
    nonempty = (e >= s).astype(np.int64)
    points_g = np.concatenate((g, g))
    points_v = np.concatenate((s, e + 1))
    delta = np.concatenate((nonempty, -nonempty))
    order, ug, uv, uid = _unique_points(points_g, points_v)
    d = np.zeros(len(uv), dtype=np.int64)
    np.add.at(d, uid, delta[order])
    cover = np.cumsum(d)  # each group's deltas sum to zero, so this stays per-group
    keep = (ug[:-1] == ug[1:]) & (cover[:-1] > 0)
    idx = np.flatnonzero(keep)
    return uv[idx], uv[idx + 1] - 1, ug[idx]


def gaps_kernel(
    s: IntArray,
    e: IntArray,
    g: IntArray,
    groups: IntArray,
    lo: IntArray,
    hi: IntArray,
) -> Triple:
    """Per group in `groups` (sorted unique), the uncovered positions within [lo, hi]."""
    rs, re, rg, _ = reduce_kernel(s, e, g, 1, True)
    gi = np.searchsorted(groups, rg)
    gi = np.minimum(gi, max(len(groups) - 1, 0))
    known = (len(groups) > 0) & (groups[gi] == rg) if len(groups) else np.zeros(len(rg), dtype=bool)
    rs, re, rg, gi = rs[known], re[known], rg[known], gi[known]
    rs, re = np.maximum(rs, lo[gi]), np.minimum(re, hi[gi])
    inside = rs <= re
    rs, re, rg = rs[inside], re[inside], rg[inside]
    gap_s_g = np.concatenate((groups, rg))
    gap_s = np.concatenate((lo, re + 1))
    gap_e_g = np.concatenate((groups, rg))
    gap_e = np.concatenate((hi, rs - 1))
    so = np.lexsort((gap_s, gap_s_g))
    eo = np.lexsort((gap_e, gap_e_g))
    out_s, out_e, out_g = gap_s[so], gap_e[eo], gap_s_g[so]
    ok = out_s <= out_e
    return out_s[ok], out_e[ok], out_g[ok]
