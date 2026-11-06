"""Set-operation kernels: position-set algebra and the parallel ``p*`` versions."""

from __future__ import annotations

import numpy as np

from ._inter import reduce_k

I64 = np.int64
Pair = tuple[np.ndarray, np.ndarray]


def _cover(bp: np.ndarray, start: np.ndarray, width: np.ndarray) -> np.ndarray:
    """Per-segment coverage count for segments ``[bp[i], bp[i+1] - 1]``."""
    ne = width > 0
    delta = np.zeros(len(bp), dtype=I64)
    np.add.at(delta, np.searchsorted(bp, start[ne]), 1)
    np.add.at(delta, np.searchsorted(bp, start[ne] + width[ne]), -1)
    return np.cumsum(delta)[:-1]


def union_k(xs: np.ndarray, xw: np.ndarray, ys: np.ndarray, yw: np.ndarray) -> Pair:
    s, w, _ = reduce_k(np.concatenate((xs, ys)), np.concatenate((xw, yw)), 1, drop_empty=True)
    return s, w


def _combine(xs: np.ndarray, xw: np.ndarray, ys: np.ndarray, yw: np.ndarray, keep: str) -> Pair:
    """Positions in both (``"both"``) or only in x (``"x_only"``), as reduced ranges."""
    nonempty = np.concatenate((xw, yw)) > 0
    if not nonempty.any() or (keep == "x_only" and not (xw > 0).any()):
        return np.array([], dtype=I64), np.array([], dtype=I64)
    allstart = np.concatenate((xs, ys))[nonempty]
    allend1 = (np.concatenate((xs + xw, ys + yw)))[nonempty]
    bp = np.unique(np.concatenate((allstart, allend1)))
    cx = _cover(bp, xs, xw) > 0
    cy = _cover(bp, ys, yw) > 0
    sel = (cx & cy) if keep == "both" else (cx & ~cy)
    idx = np.flatnonzero(sel)
    s = bp[idx]
    w = bp[idx + 1] - bp[idx]
    s, w, _ = reduce_k(s, w, 1, drop_empty=True)  # merge pieces that touch
    return s, w


def intersect_k(xs: np.ndarray, xw: np.ndarray, ys: np.ndarray, yw: np.ndarray) -> Pair:
    return _combine(xs, xw, ys, yw, "both")


def setdiff_k(xs: np.ndarray, xw: np.ndarray, ys: np.ndarray, yw: np.ndarray) -> Pair:
    return _combine(xs, xw, ys, yw, "x_only")


def _same_length(xs: np.ndarray, ys: np.ndarray) -> None:
    if len(xs) != len(ys):
        raise ValueError("'x' and 'y' must have the same length")


def punion_k(xs: np.ndarray, xw: np.ndarray, ys: np.ndarray, yw: np.ndarray, fill_gap: bool) -> Pair:
    _same_length(xs, ys)
    xe, ye = xs + xw - 1, ys + yw - 1
    if not fill_gap and len(xs):
        gap = np.maximum(xs, ys) - np.minimum(xe, ye) - 1
        if np.any(gap > 0):
            raise ValueError(
                "some pair of ranges have a gap within the 2 members of the pair; "
                "use fill_gap=True to enforce their union by filling the gap"
            )
    s = np.minimum(xs, ys)
    e = np.maximum(xe, ye)
    return s, e - s + 1


def pintersect_k(
    xs: np.ndarray, xw: np.ndarray, ys: np.ndarray, yw: np.ndarray, resolve_empty: str
) -> Pair:
    """Elementwise intersection, mirroring R's treatment of empty results."""
    _same_length(xs, ys)
    if resolve_empty not in ("none", "max.start", "start.x"):
        raise ValueError("'resolve_empty' must be 'none', 'max.start' or 'start.x'")
    xe, ye = xs + xw - 1, ys + yw - 1
    a_start = np.maximum(xs, ys)
    a_width = np.minimum(xe, ye) - a_start + 1
    keep_x = (xw == 0) & (((xs >= ys) & (xs <= ye)) | ((xs == ys) & (yw == 0)))
    a_start = np.where(keep_x, xs, a_start)
    a_width = np.where(keep_x, 0, a_width)
    keep_y = (yw == 0) & (ys >= xs) & (ys <= xe)
    a_start = np.where(keep_y, ys, a_start)
    a_width = np.where(keep_y, 0, a_width)
    check = (a_width < 0) & ~(keep_x | keep_y)
    if check.any():
        if resolve_empty == "none":
            raise ValueError(
                "some intersections produce ambiguous empty ranges; use 'resolve_empty' to resolve them"
            )
        a_width = np.where(check, 0, a_width)
        if resolve_empty == "start.x":
            a_start = np.where(check, xs, a_start)
    return a_start.astype(I64), a_width.astype(I64)


def psetdiff_k(xs: np.ndarray, xw: np.ndarray, ys: np.ndarray, yw: np.ndarray) -> Pair:
    _same_length(xs, ys)
    xe, ye = xs + xw - 1, ys + yw - 1
    if np.any((ys > xs) & (ye < xe)):
        raise ValueError(
            "some ranges in 'y' have their end points strictly inside the range in 'x' "
            "that they need to be subtracted from; cannot subtract them"
        )
    start2 = np.maximum(xs, ys)
    end2 = np.minimum(xe, ye)
    ii = start2 <= end2
    jj = end2 == xe
    new_end = np.where(ii & jj, start2 - 1, xe)
    new_start = np.where(ii & ~jj, end2 + 1, xs)
    return new_start.astype(I64), (new_end - new_start + 1).astype(I64)
