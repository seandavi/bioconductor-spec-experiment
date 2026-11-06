"""Set operations on start/width arrays.

Sequence ops (union/intersect/setdiff) work on the non-empty coverage and
return normalised ranges (sorted, non-adjacent, no zero-width). Parallel ops
(p*) are element-wise on equal-length inputs.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from . import inter

Coords = tuple[np.ndarray, np.ndarray]


def _combine(x: Coords, y: Coords, keep: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> Coords:
    xs, xw, _ = inter.reduce_runs(*x, 1, True)
    ys, yw, _ = inter.reduce_runs(*y, 1, True)
    bp = np.unique(np.concatenate([xs, xs + xw, ys, ys + yw]))
    if not len(bp):
        return inter.EMPTY, inter.EMPTY

    def covered(s: np.ndarray, w: np.ndarray) -> np.ndarray:
        delta = np.zeros(len(bp) + 1, dtype=np.int64)
        np.add.at(delta, np.searchsorted(bp, s), 1)
        np.add.at(delta, np.searchsorted(bp, s + w), -1)
        return np.cumsum(delta)[:-1] > 0

    sel = np.flatnonzero(keep(covered(xs, xw), covered(ys, yw))[:-1])
    seg_start = bp[sel]
    seg_width = bp[sel + 1] - bp[sel]
    rs, rw, _ = inter.reduce_runs(seg_start, seg_width, 1, True)
    return rs, rw


def union(x: Coords, y: Coords) -> Coords:
    return _combine(x, y, lambda a, b: a | b)


def intersect(x: Coords, y: Coords) -> Coords:
    return _combine(x, y, lambda a, b: a & b)


def setdiff(x: Coords, y: Coords) -> Coords:
    return _combine(x, y, lambda a, b: a & ~b)


def _same_length(x: Coords, y: Coords) -> None:
    if len(x[0]) != len(y[0]):
        raise ValueError("'x' and 'y' must have the same length")


def punion(x: Coords, y: Coords, fill_gap: bool) -> Coords:
    _same_length(x, y)
    xe, ye = x[0] + x[1] - 1, y[0] + y[1] - 1
    lo, hi = np.minimum(x[0], y[0]), np.maximum(xe, ye)
    if not fill_gap and (np.maximum(x[0], y[0]) > np.minimum(xe, ye) + 1).any():
        raise ValueError("some pairs have a gap; use fill_gap=True")
    return lo, hi - lo + 1


def pintersect(x: Coords, y: Coords, resolve_empty: str) -> Coords:
    if resolve_empty not in ("none", "start", "max.start"):
        raise ValueError("'resolve_empty' must be 'none', 'start' or 'max.start'")
    _same_length(x, y)
    xe, ye = x[0] + x[1] - 1, y[0] + y[1] - 1
    s, e = np.maximum(x[0], y[0]), np.minimum(xe, ye)
    ambiguous = e < s - 1
    if ambiguous.any() and resolve_empty == "none":
        raise ValueError("some intersections are ambiguous empty ranges; set resolve_empty")
    if resolve_empty == "start":
        s = np.where(ambiguous, x[0], s)
        e = np.where(ambiguous, x[0] - 1, e)
    elif resolve_empty == "max.start":
        e = np.where(ambiguous, s - 1, e)
    return s, e - s + 1


def psetdiff(x: Coords, y: Coords) -> Coords:
    _same_length(x, y)
    xs, xe = x[0], x[0] + x[1] - 1
    ys, ye = y[0], y[0] + y[1] - 1
    if ((ys > xs) & (ye < xe)).any():
        raise ValueError("some ranges in 'y' lie strictly inside 'x'; cannot subtract")
    disjoint = (ye < xs) | (ys > xe)
    covers_all = (ys <= xs) & (ye >= xe)
    new_s = np.where(~disjoint & (ys <= xs), ye + 1, xs)
    new_e = np.where(~disjoint & (ye >= xe), ys - 1, xe)
    new_s = np.where(~disjoint & covers_all, xs, new_s)
    new_e = np.where(~disjoint & covers_all, xs - 1, new_e)
    return new_s, new_e - new_s + 1
