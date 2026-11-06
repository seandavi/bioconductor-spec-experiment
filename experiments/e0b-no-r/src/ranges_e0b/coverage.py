"""Coverage kernel."""

from typing import Any

import numpy as np

from ._sew import IntArray
from .rle import Rle


def coverage_kernel(
    start: IntArray, end: IntArray, shift: int = 0, width: int | None = None, weight: Any = 1
) -> Rle:
    """Coverage over positions 1..width (default: max shifted end, 0 if none).

    Empty ranges are ignored; positions outside [1, width] are clipped.
    """
    w = np.asarray(weight)
    if w.ndim == 0:
        w = np.full(len(start), w)
    if len(w) != len(start):
        raise ValueError("weight must be scalar or one per range")
    if w.dtype.kind not in "iuf":
        raise ValueError("weight must be numeric")
    dtype = np.float64 if w.dtype.kind == "f" else np.int64
    w = w.astype(dtype, copy=False)
    s, e = start + shift, end + shift
    nonempty = e >= s
    s, e, w = s[nonempty], e[nonempty], w[nonempty]
    if width is None:
        width = max(int(e.max()), 0) if len(e) else 0
    if width < 0:
        raise ValueError("width must be non-negative")
    s, e = np.maximum(s, 1), np.minimum(e, width)
    keep = s <= e
    s, e, w = s[keep], e[keep], w[keep]
    diff = np.zeros(width + 1, dtype=dtype)
    np.add.at(diff, s - 1, w)
    np.add.at(diff, e, -w)
    return Rle.from_array(np.cumsum(diff)[:width])
