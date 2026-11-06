"""findOverlaps kernel on (start, width) arrays. Returns all hit pairs, sorted by (query, subject)."""

from __future__ import annotations

import numpy as np

from ._util import IntArray

TYPES = ("any", "start", "end", "within", "equal")


def check_overlap_args(maxgap: int, minoverlap: int, type: str) -> None:
    if type not in TYPES:
        raise ValueError(f"type must be one of {TYPES}")
    if maxgap < -1:
        raise ValueError("maxgap must be >= -1")
    if minoverlap < 0:
        raise ValueError("minoverlap cannot be negative")
    if type == "any" and maxgap != -1 and minoverlap != 0:
        raise ValueError("when type is 'any', at least one of maxgap and minoverlap must be at its default")


def overlap_pairs(
    qs: IntArray,
    qw: IntArray,
    ss: IntArray,
    sw: IntArray,
    maxgap: int = -1,
    minoverlap: int = 0,
    type: str = "any",
) -> tuple[IntArray, IntArray]:
    """All (query, subject) index pairs satisfying the R overlap rule for `type`.

    any:    gap(q, s) <= maxgap (gap = -1 when overlapping; adjacent = 0)
    start:  |start(q) - start(s)| <= max(maxgap, 0); no overlap needed
    end:    |end(q) - end(s)| <= max(maxgap, 0)
    equal:  both of the above
    within: q inside s; if maxgap > 0 also width(s) - width(q) <= maxgap
    Every type additionally needs `minoverlap` shared positions.
    """
    check_overlap_args(maxgap, minoverlap, type)
    nq, ns = len(qs), len(ss)
    empty = np.zeros(0, dtype=np.int64)
    if nq == 0 or ns == 0:
        return empty, empty
    qe = qs + qw - 1
    se = ss + sw - 1
    by_start = np.argsort(ss, kind="stable")
    s_start = ss[by_start]
    s_end = se[by_start]
    prefix_end = np.maximum.accumulate(s_end)
    by_end = np.argsort(se, kind="stable")
    e_sorted = se[by_end]
    m = max(maxgap, 0)
    q_out: list[IntArray] = []
    s_out: list[IntArray] = []
    for i in range(nq):
        if type == "any":
            hi = np.searchsorted(s_start, qe[i] + maxgap + 1, side="right")
            lo = np.searchsorted(prefix_end, qs[i] - maxgap - 1, side="left")
            cand = by_start[lo:hi]
            cand = cand[se[cand] >= qs[i] - maxgap - 1]
        elif type in ("start", "equal"):
            lo = np.searchsorted(s_start, qs[i] - m, side="left")
            hi = np.searchsorted(s_start, qs[i] + m, side="right")
            cand = by_start[lo:hi]
            if type == "equal":
                cand = cand[np.abs(se[cand] - qe[i]) <= m]
        elif type == "end":
            lo = np.searchsorted(e_sorted, qe[i] - m, side="left")
            hi = np.searchsorted(e_sorted, qe[i] + m, side="right")
            cand = by_end[lo:hi]
        else:  # within
            hi = np.searchsorted(s_start, qs[i], side="right")
            lo = np.searchsorted(prefix_end, qe[i], side="left")
            cand = by_start[lo:hi]
            cand = cand[(se[cand] >= qe[i]) & (ss[cand] <= qs[i])]
            if maxgap > 0:
                cand = cand[sw[cand] - qw[i] <= maxgap]
        if minoverlap > 0 and len(cand):
            shared = np.minimum(qe[i], se[cand]) - np.maximum(qs[i], ss[cand]) + 1
            cand = cand[shared >= minoverlap]
        if len(cand):
            cand = np.sort(cand)
            q_out.append(np.full(len(cand), i, dtype=np.int64))
            s_out.append(cand.astype(np.int64))
    if not q_out:
        return empty, empty
    return np.concatenate(q_out), np.concatenate(s_out)
