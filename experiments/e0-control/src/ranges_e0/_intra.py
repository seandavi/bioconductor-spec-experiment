"""Intra-range kernels: pure functions on ``(start, width)`` int64 arrays.

Each returns a new ``(start, width)`` pair of the same length.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ._validate import as_int_array, recycle

Pair = tuple[np.ndarray, np.ndarray]


def _opt(x: Any, n: int, name: str) -> np.ndarray | None:
    return None if x is None else recycle(as_int_array(x, name), n, name)


def shift_k(start: np.ndarray, width: np.ndarray, shift: Any) -> Pair:
    sh = recycle(as_int_array(shift, "shift"), len(start), "shift") if len(start) else np.array([], dtype=np.int64)
    return start + sh, width.copy()


def narrow_k(start: np.ndarray, width: np.ndarray, s: Any, e: Any, w: Any) -> Pair:
    """R's ``narrow``: ``s``/``e`` are 1-based offsets *within each range*;
    negative values count from the end of the range. Without
    ``allow.nonnarrowing`` the result must lie inside the original range."""
    n = len(start)
    if n == 0:
        return start.copy(), width.copy()
    sa, ea, wa = _opt(s, n, "start"), _opt(e, n, "end"), _opt(w, n, "width")
    if wa is not None and sa is not None and ea is not None:
        raise ValueError("either the supplied start or the supplied end (but not both) must be NA when width is supplied")
    if wa is not None and sa is None and ea is None:
        raise ValueError("either the supplied start or the supplied end (but not both) must be NA when width is supplied")
    ref = width
    if sa is not None:
        sa = np.where(sa < 0, ref + sa + 1, sa)
    if ea is not None:
        ea = np.where(ea < 0, ref + ea + 1, ea)
    if wa is not None:
        if np.any(wa < 0):
            raise ValueError("'width' must be non-negative")
        if sa is not None:
            ea = sa + wa - 1
        else:
            assert ea is not None
            sa = ea - wa + 1
    if sa is None:
        sa = np.ones(n, dtype=np.int64)
    if ea is None:
        ea = ref.copy()
    if np.any(sa < 1):
        raise ValueError("'allow.nonnarrowing' is FALSE and the supplied start is < 1")
    if np.any(sa > ref + 1):
        raise ValueError("'allow.nonnarrowing' is FALSE and the supplied start is > refwidth + 1")
    if np.any(ea > ref):
        raise ValueError("'allow.nonnarrowing' is FALSE and the supplied end is > refwidth")
    if np.any(ea < sa - 1):
        raise ValueError("the supplied end is < start - 1")
    return start + sa - 1, ea - sa + 1


def resize_k(start: np.ndarray, width: np.ndarray, new_width: Any, fix: np.ndarray) -> Pair:
    """``fix`` is an object array of ``"start"`` / ``"end"`` / ``"center"``."""
    n = len(start)
    if n == 0:
        return start.copy(), width.copy()
    nw = recycle(as_int_array(new_width, "width"), n, "width")
    if np.any(nw < 0):
        raise ValueError("'width' values must be non-negative")
    end = start + width - 1
    new_start = start.copy()
    is_end = fix == "end"
    is_ctr = fix == "center"
    new_start = np.where(is_end, end - nw + 1, new_start)
    new_start = np.where(is_ctr, start + (width - nw) // 2, new_start)
    return new_start, nw.copy()


def as_fix(fix: Any, n: int) -> np.ndarray:
    """Validate/recycle the ``fix`` argument to an object array of length ``n``."""
    arr = np.empty(n, dtype=object)
    vals = [fix] if isinstance(fix, str) else list(fix)
    if len(vals) not in (1, n) and n:
        raise ValueError("'fix' must have length 1 or the length of the object")
    for v in vals:
        if v not in ("start", "end", "center"):
            raise ValueError("'fix' must be 'start', 'end' or 'center'")
    if n:
        arr[:] = vals * n if len(vals) == 1 else vals
    return arr


def flank_k(start: np.ndarray, width: np.ndarray, flank_width: Any, at_start: np.ndarray, both: bool) -> Pair:
    """``at_start`` is a bool array: flank the start (True) or the end (False)."""
    n = len(start)
    if n == 0:
        return start.copy(), width.copy()
    fw = recycle(as_int_array(flank_width, "width"), n, "width")
    end = start + width - 1
    if both:
        fw = np.abs(fw)
        return np.where(at_start, start - fw, end - fw + 1), 2 * fw
    new_start = np.where(
        at_start,
        np.where(fw >= 0, start - fw, start),
        np.where(fw >= 0, end + 1, end + fw + 1),
    )
    return new_start, np.abs(fw)


def promoters_k(start: np.ndarray, width: np.ndarray, upstream: Any, downstream: Any, neg: np.ndarray | None = None) -> Pair:
    """Window around the TSS: ``[tss - upstream, tss + downstream - 1]``.

    ``neg`` marks ranges whose TSS is their *end* (minus strand): there the
    window is ``[end - downstream + 1, end + upstream]``."""
    n = len(start)
    if n == 0:
        return start.copy(), width.copy()
    up = recycle(as_int_array(upstream, "upstream"), n, "upstream")
    dn = recycle(as_int_array(downstream, "downstream"), n, "downstream")
    if up.min() < 0 or dn.min() < 0:
        raise ValueError("'upstream' and 'downstream' must be integers >= 0")
    end = start + width - 1
    if neg is None:
        neg = np.zeros(n, dtype=bool)
    new_start = np.where(neg, end - dn + 1, start - up)
    return new_start, up + dn
