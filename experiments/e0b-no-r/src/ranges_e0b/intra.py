"""Intra-range kernels on (start, width) arrays. Each returns (start, width)."""

from typing import Any

import numpy as np

from ._sew import BoolArray, IntArray, per_range

FIX_START, FIX_END, FIX_CENTER = 0, 1, 2
_FIX_CODES = {"start": FIX_START, "end": FIX_END, "center": FIX_CENTER}


def parse_fix(fix: str | list[str] | tuple[str, ...] | np.ndarray, n: int) -> IntArray:
    """`fix` as scalar or per-range sequence -> int codes of length n."""
    vals = [fix] if isinstance(fix, str) else list(fix)
    try:
        codes = np.array([_FIX_CODES[v] for v in vals], dtype=np.int64)
    except KeyError as err:
        raise ValueError(f"fix must be 'start', 'end' or 'center', got {err.args[0]!r}") from None
    if len(codes) == 1:
        return np.repeat(codes, n)
    if len(codes) != n:
        raise ValueError(f"fix must have length 1 or {n}")
    return codes


def shift(start: IntArray, width: IntArray, k: Any) -> tuple[IntArray, IntArray]:
    return start + per_range(k, len(start), "shift"), width


def narrow(
    start: IntArray, width: IntArray, a_start: Any, a_end: Any, a_width: Any
) -> tuple[IntArray, IntArray]:
    """Relative start/end/width; negative start/end count back from range end."""
    n = len(start)
    if a_start is not None and a_end is not None and a_width is not None:
        raise ValueError("at most two of start, end, width may be given")
    rs: IntArray | None = None
    re: IntArray | None = None
    if a_start is not None:
        v = per_range(a_start, n, "start")
        rs = np.where(v < 0, width + v + 1, v)
    if a_end is not None:
        v = per_range(a_end, n, "end")
        re = np.where(v < 0, width + v + 1, v)
    rw = None if a_width is None else per_range(a_width, n, "width")
    if rw is not None and np.any(rw < 0):
        raise ValueError("width must be non-negative")
    if rs is None and re is None and rw is None:
        return start, width
    if rs is not None and re is not None:
        pass
    elif rs is not None:
        re = rs + rw - 1 if rw is not None else width.copy()
    elif re is not None:
        rs = re - rw + 1 if rw is not None else np.ones(n, dtype=np.int64)
    else:
        assert rw is not None
        rs = np.ones(n, dtype=np.int64)  # width alone keeps the leading positions
        re = rw.copy()
    assert rs is not None and re is not None
    if np.any(rs < 1) or np.any(re > width) or np.any(rs > re + 1):
        raise ValueError("narrowed range falls outside the original range")
    return start + rs - 1, re - rs + 1


def resize(start: IntArray, width: IntArray, new_width: Any, fix: IntArray) -> tuple[IntArray, IntArray]:
    n = len(start)
    w2 = per_range(new_width, n, "width")
    if np.any(w2 < 0):
        raise ValueError("width must be non-negative")
    end = start + width - 1
    s2 = np.where(
        fix == FIX_START,
        start,
        np.where(fix == FIX_END, end - w2 + 1, start + (width - w2) // 2),
    )
    return s2, w2


def flank(
    start: IntArray, width: IntArray, flank_width: Any, at_start: BoolArray, both: bool
) -> tuple[IntArray, IntArray]:
    """Flanks at range start (`at_start`) or end; negative width flanks inward."""
    n = len(start)
    w = per_range(flank_width, n, "width")
    end = start + width - 1
    aw = np.abs(w)
    if both:
        s2 = np.where(at_start, start - aw, end - aw + 1)
        return s2, 2 * aw
    s2 = np.where(
        at_start,
        np.where(w >= 0, start - w, start),
        np.where(w >= 0, end + 1, end - aw + 1),
    )
    return s2, aw


def promoters(
    start: IntArray, width: IntArray, upstream: Any, downstream: Any, minus: BoolArray
) -> tuple[IntArray, IntArray]:
    n = len(start)
    up = per_range(upstream, n, "upstream")
    down = per_range(downstream, n, "downstream")
    if np.any(up < 0) or np.any(down < 0):
        raise ValueError("upstream and downstream must be non-negative")
    end = start + width - 1
    s2 = np.where(minus, end - down + 1, start - up)
    return s2, up + down


def as_bool_array(x: Any, n: int, name: str) -> BoolArray:
    a = np.asarray(x)
    if a.ndim == 0:
        a = a.reshape(1)
    if a.dtype.kind != "b":
        raise ValueError(f"{name} must be boolean")
    if len(a) == 1:
        return np.repeat(a, n)
    if len(a) != n:
        raise ValueError(f"{name} must have length 1 or {n}")
    return a

