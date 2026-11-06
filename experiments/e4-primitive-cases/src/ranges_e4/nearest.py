"""nearest / precede / follow / distance kernels over start/width arrays.

distance(x, y) = max(0, max(xs, ys) - min(xe, ye) - 1): 0 for overlapping or
adjacent ranges. ``nearest`` minimises it over all subjects, so overlapping
and adjacent subjects tie. ``precede`` takes subjects with start > end(x) and
the smallest start; ``follow`` takes subjects with end < start(x) and the
largest end. Ties are all returned, sorted by (query, subject); callers pick
first/last/all. Zero-width ranges use their ``[s, s-1]`` coordinates.
"""

from __future__ import annotations

import numpy as np

Pairs = tuple[np.ndarray, np.ndarray]
EMPTY_PAIRS: Pairs = (np.zeros(0, np.int64), np.zeros(0, np.int64))


def distance(xs: np.ndarray, xw: np.ndarray, ys: np.ndarray, yw: np.ndarray) -> np.ndarray:
    xe, ye = xs + xw - 1, ys + yw - 1
    return np.maximum(0, np.maximum(xs, ys) - np.minimum(xe, ye) - 1)


def _collect(per_query: list[np.ndarray]) -> Pairs:
    q = [np.full(len(h), i, dtype=np.int64) for i, h in enumerate(per_query) if len(h)]
    if not q:
        return EMPTY_PAIRS
    return np.concatenate(q), np.concatenate([h for h in per_query if len(h)]).astype(np.int64)


def nearest_pairs(qs: np.ndarray, qw: np.ndarray, ss: np.ndarray, sw: np.ndarray) -> Pairs:
    if not len(ss):
        return EMPTY_PAIRS
    out = []
    for i in range(len(qs)):
        d = distance(qs[i], qw[i], ss, sw)
        out.append(np.flatnonzero(d == d.min()))
    return _collect(out)


def precede_pairs(qs: np.ndarray, qw: np.ndarray, ss: np.ndarray, sw: np.ndarray) -> Pairs:
    qe = qs + qw - 1
    out = []
    for i in range(len(qs)):
        cand = ss > qe[i]
        out.append(np.flatnonzero(cand & (ss == ss[cand].min())) if cand.any() else np.zeros(0, np.int64))
    return _collect(out)


def follow_pairs(qs: np.ndarray, qw: np.ndarray, ss: np.ndarray, sw: np.ndarray) -> Pairs:
    se = ss + sw - 1
    out = []
    for i in range(len(qs)):
        cand = se < qs[i]
        out.append(np.flatnonzero(cand & (se == se[cand].max())) if cand.any() else np.zeros(0, np.int64))
    return _collect(out)


def ties_at_minimum(dist: np.ndarray, eligible: np.ndarray) -> np.ndarray:
    """Positions of eligible entries with the smallest distance (all ties, ascending)."""
    if not eligible.any():
        return np.zeros(0, np.int64)
    d = np.where(eligible, dist, np.iinfo(np.int64).max)
    return np.flatnonzero(d == d.min())


def collect(per_query: list[np.ndarray]) -> Pairs:
    return _collect(per_query)
