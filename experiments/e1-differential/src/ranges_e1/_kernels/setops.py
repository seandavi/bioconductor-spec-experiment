"""Set operations: normalising (union / intersect / setdiff) and parallel (p*) versions."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from . import inter

_SHIFT = 1 << 33
_OFF = 1 << 32
_Triple = tuple[np.ndarray, np.ndarray, np.ndarray]


def _lift(group: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Map (group, position) onto one axis so groups never touch."""
    return group * _SHIFT + pos + _OFF


def _covered(starts: np.ndarray, ends1: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """pts[i] inside one of the sorted, disjoint [starts, ends1) ranges."""
    idx = np.searchsorted(starts, pts, side="right") - 1
    ok = idx >= 0
    ok[ok] = pts[ok] < ends1[idx[ok]]
    return ok


def combine(
    gx: np.ndarray, sx: np.ndarray, wx: np.ndarray,
    gy: np.ndarray, sy: np.ndarray, wy: np.ndarray,
    op: Callable[[np.ndarray, np.ndarray], np.ndarray],
) -> _Triple:
    """Per-group set operation on the positions covered by x and y.

    Empty ranges cover nothing and are ignored. `op(in_x, in_y)` picks segments.
    Returns (group, start, width) with adjacent pieces merged.
    """
    nx, ny = wx > 0, wy > 0
    rx = inter.reduce(gx[nx], sx[nx], wx[nx], 1)
    ry = inter.reduce(gy[ny], sy[ny], wy[ny], 1)
    lx, ly = _lift(rx.group, rx.start), _lift(ry.group, ry.start)
    ex, ey = lx + rx.width, ly + ry.width
    bps = np.unique(np.concatenate([lx, ex, ly, ey]))
    empty = np.empty(0, dtype=np.int64)
    if len(bps) < 2:
        return empty, empty, empty
    seg = bps[:-1]
    sel = op(_covered(lx, ex, seg), _covered(ly, ey, seg))
    a = seg[sel]
    b = bps[1:][sel]
    if len(a) == 0:
        return empty, empty, empty
    new = np.ones(len(a), dtype=bool)
    new[1:] = a[1:] != b[:-1]
    firsts = np.flatnonzero(new)
    lasts = np.append(firsts[1:], len(a)) - 1
    ls, le = a[firsts], b[lasts]
    return (ls >> 33), (ls & (_SHIFT - 1)) - _OFF, le - ls




def union(gx: np.ndarray, sx: np.ndarray, wx: np.ndarray, gy: np.ndarray, sy: np.ndarray, wy: np.ndarray) -> _Triple:
    return combine(gx, sx, wx, gy, sy, wy, lambda a, b: a | b)


def intersect(gx: np.ndarray, sx: np.ndarray, wx: np.ndarray, gy: np.ndarray, sy: np.ndarray, wy: np.ndarray) -> _Triple:
    return combine(gx, sx, wx, gy, sy, wy, lambda a, b: a & b)


def setdiff(gx: np.ndarray, sx: np.ndarray, wx: np.ndarray, gy: np.ndarray, sy: np.ndarray, wy: np.ndarray) -> _Triple:
    return combine(gx, sx, wx, gy, sy, wy, lambda a, b: a & ~b)


# ---- parallel versions ------------------------------------------------------


def punion(
    sx: np.ndarray, wx: np.ndarray, sy: np.ndarray, wy: np.ndarray, fill_gap: bool
) -> tuple[np.ndarray, np.ndarray]:
    ex, ey = sx + wx - 1, sy + wy - 1
    gap = np.maximum(sx, sy) - np.minimum(ex, ey) - 1
    if not fill_gap and (gap > 0).any():
        raise ValueError(
            "some pair of ranges have a gap within the 2 members of the pair. "
            "Use 'fill.gap=TRUE' to enforce their union by filling the gap."
        )
    s = np.minimum(sx, sy)
    return s, np.maximum(ex, ey) - s + 1


def pintersect(
    sx: np.ndarray, wx: np.ndarray, sy: np.ndarray, wy: np.ndarray, resolve_empty: str
) -> tuple[np.ndarray, np.ndarray]:
    ex, ey = sx + wx - 1, sy + wy - 1
    s = np.maximum(sx, sy)
    w = np.minimum(ex, ey) - s + 1
    amb = w < 0
    if amb.any():
        if resolve_empty == "none":
            raise ValueError(
                "some intersections produce ambiguous empty ranges. "
                "Use argument 'resolve.empty' to resolve them."
            )
        if resolve_empty not in ("max.start", "start.x"):
            raise ValueError("'arg' should be one of \"none\", \"max.start\", \"start.x\"")
        s = np.where(amb & (resolve_empty == "start.x"), sx, s)
        w = np.where(amb, 0, w)
    return s, w


def psetdiff(
    sx: np.ndarray, wx: np.ndarray, sy: np.ndarray, wy: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    ex, ey = sx + wx - 1, sy + wy - 1
    if ((sy > sx) & (ey < ex)).any():
        raise ValueError(
            "some ranges in 'y' have their end points strictly inside the range in 'x' "
            "that they need to be subtracted from. Cannot subtract them."
        )
    left_w = np.maximum(np.minimum(ex, sy - 1) - sx + 1, 0)
    rs = np.maximum(sx, ey + 1)
    rw = ex - rs + 1
    right_ok = rw > 0
    start = np.where(sy > sx, sx, np.where(right_ok, rs, sx))
    width = np.where(sy > sx, left_w, np.where(right_ok, rw, 0))
    return start, width
