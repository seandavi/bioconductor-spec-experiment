"""Set-operation kernels (grouped) and parallel (pairwise) kernels."""

import numpy as np

from ._sew import IntArray
from .inter import Triple, _unique_points, reduce_kernel


def union_kernel(sx: IntArray, ex: IntArray, gx: IntArray, sy: IntArray, ey: IntArray, gy: IntArray) -> Triple:
    s, e, g, _ = reduce_kernel(
        np.concatenate((sx, sy)), np.concatenate((ex, ey)), np.concatenate((gx, gy)), 1, True
    )
    return s, e, g


def sweep_kernel(
    sx: IntArray, ex: IntArray, gx: IntArray, sy: IntArray, ey: IntArray, gy: IntArray, op: str
) -> Triple:
    """Positions covered by x and (op="intersect") or not (op="setdiff") by y, per group.

    Same result as R's gaps()-based definitions: normalised (reduced, non-empty) output.
    """
    if op not in ("intersect", "setdiff"):
        raise ValueError(op)
    rxs, rxe, rxg, _ = reduce_kernel(sx, ex, gx, 1, True)
    rys, rye, ryg, _ = reduce_kernel(sy, ey, gy, 1, True)
    nx, ny = len(rxs), len(rys)
    if nx == 0:
        return rxs, rxe, rxg
    g = np.concatenate((rxg, rxg, ryg, ryg))
    v = np.concatenate((rxs, rxe + 1, rys, rye + 1))
    dx = np.concatenate((np.ones(nx, np.int64), -np.ones(nx, np.int64), np.zeros(2 * ny, np.int64)))
    dy = np.concatenate((np.zeros(2 * nx, np.int64), np.ones(ny, np.int64), -np.ones(ny, np.int64)))
    order, ug, uv, uid = _unique_points(g, v)
    ax = np.zeros(len(uv), np.int64)
    ay = np.zeros(len(uv), np.int64)
    np.add.at(ax, uid, dx[order])
    np.add.at(ay, uid, dy[order])
    cx, cy = np.cumsum(ax), np.cumsum(ay)
    keep = (cx > 0) & ((cy > 0) if op == "intersect" else (cy == 0))
    idx = np.flatnonzero(keep)
    s, e, gg, _ = reduce_kernel(uv[idx], uv[idx + 1] - 1, ug[idx], 1, True)
    return s, e, gg


# --- parallel ----------------------------------------------------------------


def _same_length(*arrays: IntArray) -> None:
    if len({len(a) for a in arrays}) != 1:
        raise ValueError("parallel operations require equal lengths")


def punion(sx: IntArray, ex: IntArray, sy: IntArray, ey: IntArray, fill_gap: bool) -> tuple[IntArray, IntArray]:
    _same_length(sx, sy)
    if not fill_gap and np.any(np.maximum(sx, sy) - np.minimum(ex, ey) - 1 > 0):
        raise ValueError("ranges are separated by a gap; use fill_gap=True")
    return np.minimum(sx, sy), np.maximum(ex, ey)


def pintersect(
    sx: IntArray, ex: IntArray, sy: IntArray, ey: IntArray, resolve_empty: str
) -> tuple[IntArray, IntArray]:
    _same_length(sx, sy)
    if resolve_empty not in ("none", "max.start", "start.x"):
        raise ValueError(f"unknown resolve_empty {resolve_empty!r}")
    s, e = np.maximum(sx, sy), np.minimum(ex, ey)
    bad = e < s - 1
    if np.any(bad):
        if resolve_empty == "none":
            raise ValueError("some ranges do not overlap; set resolve_empty")
        s_fix = s if resolve_empty == "max.start" else sx
        s = np.where(bad, s_fix, s)
        e = np.where(bad, s_fix - 1, e)
    return s, e


def psetdiff(sx: IntArray, ex: IntArray, sy: IntArray, ey: IntArray) -> tuple[IntArray, IntArray]:
    """x minus y per pair; y strictly inside x would split it and is an error."""
    _same_length(sx, sy)
    hit = (sy <= ex) & (ey >= sx) & (ey >= sy)
    covers = hit & (sy <= sx) & (ey >= ex)
    cut_left = hit & (sy <= sx) & ~covers
    cut_right = hit & (ey >= ex) & ~covers & ~cut_left
    split = hit & ~covers & ~cut_left & ~cut_right
    if np.any(split):
        raise ValueError("y lies strictly inside x; result would need two ranges")
    s = np.where(cut_left, ey + 1, sx)
    e = np.where(cut_right, sy - 1, ex)
    s = np.where(covers, sx, s)
    e = np.where(covers, sx - 1, e)
    return s, e
