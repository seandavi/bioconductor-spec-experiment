"""precede / follow / nearest / distance.

One driver serves IRanges and GRanges: callers pass per-range seqlevel codes
and strand codes (0 `+`, 1 `-`, 2 `*`); IRanges passes all-0 seq and all-2 strand.
"""

import numpy as np

from ._sew import NO_HIT, IntArray, empty_int
from .overlaps import overlap_pairs

KINDS = ("precede", "follow", "nearest")
MINUS, STAR = 1, 2


def distance(sx: IntArray, ex: IntArray, sy: IntArray, ey: IntArray) -> IntArray:
    """Parallel distance; overlapping or adjacent ranges are 0."""
    if len(sx) != len(sy):
        raise ValueError("distance requires equal lengths")
    return np.maximum(np.maximum(sx, sy) - np.minimum(ex, ey) - 1, 0)


def precede_pairs(qe: IntArray, ss: IntArray, exclude: IntArray | None = None) -> tuple[IntArray, IntArray]:
    """For each query, all subjects with the smallest start > query end.

    `exclude[q]` is a subject index never returned for query q (self-matches).
    """
    nq, ns = len(qe), len(ss)
    if nq == 0 or ns == 0:
        return empty_int(), empty_int()
    order = np.lexsort((np.arange(ns), ss))
    sorted_s = ss[order]
    pos = np.searchsorted(sorted_s, qe, side="right")
    pending = np.arange(nq)
    out_q: list[IntArray] = []
    out_s: list[IntArray] = []
    while len(pending):
        p = pos[pending]
        valid = p < ns
        pending, p = pending[valid], p[valid]
        if not len(pending):
            break
        gend = np.searchsorted(sorted_s, sorted_s[p], side="right")
        c = gend - p
        q_idx = np.repeat(pending, c)
        offs = np.arange(int(c.sum())) - np.repeat(np.cumsum(c) - c, c)
        s_idx = order[np.repeat(p, c) + offs]
        if exclude is None:
            out_q.append(q_idx)
            out_s.append(s_idx)
            break
        ok = s_idx != exclude[q_idx]
        out_q.append(q_idx[ok])
        out_s.append(s_idx[ok])
        found = np.zeros(nq, dtype=bool)
        found[q_idx[ok]] = True
        redo = ~found[pending]  # only self was in the group: move on to the next start
        pos[pending[redo]] = gend[redo]
        pending = pending[redo]
    if not out_q:
        return empty_int(), empty_int()
    q, s = np.concatenate(out_q), np.concatenate(out_s)
    o = np.lexsort((s, q))
    return q[o], s[o]


def nearest_pairs(
    qs: IntArray, qe: IntArray, ss: IntArray, se: IntArray, exclude: IntArray | None
) -> tuple[IntArray, IntArray]:
    """All nearest subjects per query; overlapping/adjacent ones (distance 0) win."""
    ns = len(ss)
    zq, zs = overlap_pairs(qs, qe, ss, se, "any", 0, 0)
    if exclude is not None:
        keep = zs != exclude[zq]
        zq, zs = zq[keep], zs[keep]
    pq, ps = precede_pairs(qe, ss, exclude)
    fq, fs = precede_pairs(-qs, -se, exclude)
    q = np.concatenate((zq, pq, fq))
    s = np.concatenate((zs, ps, fs))
    d = np.concatenate((np.zeros(len(zq), np.int64), ss[ps] - qe[pq] - 1, qs[fq] - se[fs] - 1))
    if len(q) == 0:
        return empty_int(), empty_int()
    best = np.full(len(qs), np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(best, q, d)
    keep = d == best[q]
    key = np.unique(q[keep] * ns + s[keep])
    return key // ns, key % ns


def nearest_family(
    kind: str,
    qs: IntArray,
    qe: IntArray,
    qseq: IntArray,
    qstrand: IntArray,
    ss: IntArray,
    se: IntArray,
    sseq: IntArray,
    sstrand: IntArray,
    ignore_strand: bool,
    exclude_self: bool,
) -> tuple[IntArray, IntArray]:
    """All (query, subject) pairs for precede / follow / nearest, 0-based, sorted."""
    if kind not in KINDS:
        raise ValueError(kind)
    out_q: list[IntArray] = []
    out_s: list[IntArray] = []
    for seq in np.unique(qseq):
        in_seq_q = qseq == seq
        in_seq_s = sseq == seq
        for st in np.unique(qstrand[in_seq_q]):
            qidx = np.flatnonzero(in_seq_q & (qstrand == st))
            if ignore_strand or st == STAR:
                sidx = np.flatnonzero(in_seq_s)
            else:
                sidx = np.flatnonzero(in_seq_s & ((sstrand == st) | (sstrand == STAR)))
            if len(qidx) == 0 or len(sidx) == 0:
                continue
            exclude = np.searchsorted(sidx, qidx) if exclude_self else None
            q_s, q_e, s_s, s_e = qs[qidx], qe[qidx], ss[sidx], se[sidx]
            if kind == "nearest":
                lq, ls = nearest_pairs(q_s, q_e, s_s, s_e, exclude)
            else:
                flip = (kind == "follow") != (st == MINUS and not ignore_strand)
                if flip:
                    q_s, q_e, s_s, s_e = -q_e, -q_s, -s_e, -s_s
                lq, ls = precede_pairs(q_e, s_s, exclude)
            out_q.append(qidx[lq])
            out_s.append(sidx[ls])
    if not out_q:
        return empty_int(), empty_int()
    q, s = np.concatenate(out_q), np.concatenate(out_s)
    o = np.lexsort((s, q))
    return q[o], s[o]


def first_per_query(qi: IntArray, si: IntArray, n_query: int) -> IntArray:
    out = np.full(n_query, NO_HIT, dtype=np.int64)
    out[qi[::-1]] = si[::-1]
    return out
