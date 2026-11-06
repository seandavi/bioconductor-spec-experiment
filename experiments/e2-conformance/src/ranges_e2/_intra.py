"""Intra-range kernels on (start, width) int64 arrays: narrow, resize, flank, promoters."""

from __future__ import annotations

import numpy as np

from ._util import IntArray, rep_len


def narrow_ranges(
    start: IntArray,
    width: IntArray,
    nstart: IntArray | None,
    nend: IntArray | None,
    nwidth: IntArray | None,
) -> tuple[IntArray, IntArray]:
    """R's narrow: coordinates are relative to each range (1-based; negatives count from the end)."""
    n = len(start)
    if sum(a is not None for a in (nstart, nend, nwidth)) == 3:
        raise ValueError("only two of start, end and width can be supplied")
    if nwidth is not None and nstart is None and nend is None:
        raise ValueError("either start or end must be supplied when width is")
    s = rep_len(nstart, n, "start") if nstart is not None else None
    e = rep_len(nend, n, "end") if nend is not None else None
    w = rep_len(nwidth, n, "width") if nwidth is not None else None
    if s is not None:
        s = np.where(s < 0, width + s + 1, s)
    if e is not None:
        e = np.where(e < 0, width + e + 1, e)
    if s is None and e is None:
        s = np.ones(n, dtype=np.int64)
        e = width.copy()
    elif s is None:
        assert e is not None
        s = e - w + 1 if w is not None else np.ones(n, dtype=np.int64)
    elif e is None:
        e = s + w - 1 if w is not None else width.copy()
    assert s is not None and e is not None
    if np.any(s < 1):
        raise ValueError("narrowing start is < 1 (allow.nonnarrowing is FALSE)")
    if np.any(e > width):
        raise ValueError("narrowing end is > width (allow.nonnarrowing is FALSE)")
    if np.any(s > width + 1):
        raise ValueError("narrowing start is > width + 1 (allow.nonnarrowing is FALSE)")
    if np.any(e - s + 1 < 0):
        raise ValueError("the supplied start/end lead to a negative width")
    return start + s - 1, e - s + 1


def resize_ranges(
    start: IntArray, width: IntArray, new_width: IntArray, fix: np.ndarray
) -> tuple[IntArray, IntArray]:
    """fix: array of 'start' | 'end' | 'center'."""
    n = len(start)
    nw = rep_len(new_width, n, "width")
    fx = rep_len(fix, n, "fix")
    if np.any(nw < 0):
        raise ValueError("width values must be non-negative")
    bad = ~np.isin(fx, ["start", "end", "center"])
    if bad.any():
        raise ValueError("fix must be one of 'start', 'end', 'center'")
    out_start = start.copy()
    is_end = fx == "end"
    out_start[is_end] = (start + width - nw)[is_end]
    is_center = fx == "center"
    out_start[is_center] = (start + np.floor_divide(width - nw, 2))[is_center]
    return out_start, nw.copy()


def flank_ranges(
    start: IntArray, width: IntArray, flank_width: IntArray, at_start: np.ndarray, both: bool
) -> tuple[IntArray, IntArray]:
    n = len(start)
    fw = rep_len(flank_width, n, "width")
    st = rep_len(at_start, n, "start")
    a = np.abs(fw)
    end_excl = start + width
    if both:
        out_start = np.where(st, start - a, end_excl - a)
        return out_start, 2 * a
    left = np.where(fw >= 0, start - a, start)  # start=TRUE
    right = np.where(fw >= 0, end_excl, end_excl - a)  # start=FALSE
    return np.where(st, left, right), a


def promoter_ranges(
    start: IntArray, width: IntArray, upstream: IntArray, downstream: IntArray, minus: np.ndarray
) -> tuple[IntArray, IntArray]:
    """`minus[i]` True mirrors the window around the range end (GRanges '-' strand)."""
    n = len(start)
    up = rep_len(upstream, n, "upstream")
    down = rep_len(downstream, n, "downstream")
    if np.any(up < 0) or np.any(down < 0):
        raise ValueError("upstream and downstream must be integers >= 0")
    end = start + width - 1
    out_start = np.where(minus, end - down + 1, start - up)
    return out_start, up + down
