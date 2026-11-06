"""precede / follow on flat arrays (R: findInterval over sorted starts / ends)."""

from __future__ import annotations

import numpy as np


def precede(qe: np.ndarray, ss: np.ndarray, select_all: bool) -> tuple[np.ndarray, np.ndarray]:
    """For each query end, the subject(s) with the smallest start strictly greater than it.

    Returns (query_idx, subject_idx); ties on start are all returned (select_all) or
    only the lowest subject index. Queries with no such subject are absent.
    """
    empty = np.empty(0, dtype=np.int64)
    if len(ss) == 0 or len(qe) == 0:
        return empty, empty
    order = np.argsort(ss, kind="stable")
    ssort = ss[order]
    i = np.searchsorted(ssort, qe, side="right")
    ok = np.flatnonzero(i < len(ss))
    pos = i[ok]
    if not select_all:
        return ok.astype(np.int64), order[pos].astype(np.int64)
    hi = np.searchsorted(ssort, ssort[pos], side="right")
    return _expand(ok, pos, hi, order)


def follow(qs: np.ndarray, se: np.ndarray, select_all: bool) -> tuple[np.ndarray, np.ndarray]:
    """For each query start, the subject(s) with the largest end strictly smaller than it.

    Ties on end: all (select_all) or the highest subject index ("last").
    """
    empty = np.empty(0, dtype=np.int64)
    if len(se) == 0 or len(qs) == 0:
        return empty, empty
    order = np.argsort(se, kind="stable")
    esort = se[order]
    i = np.searchsorted(esort, qs - 1, side="right")
    ok = np.flatnonzero(i > 0)
    last = i[ok] - 1
    if not select_all:
        return ok.astype(np.int64), order[last].astype(np.int64)
    lo = np.searchsorted(esort, esort[last], side="left")
    return _expand(ok, lo, last + 1, order)


def _expand(q: np.ndarray, lo: np.ndarray, hi: np.ndarray, order: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    counts = hi - lo
    qi = np.repeat(q, counts).astype(np.int64)
    offs = np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
    return qi, order[np.repeat(lo, counts) + offs].astype(np.int64)
