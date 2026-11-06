"""Coverage kernel: sparse event sweep producing run-length encoded depth."""

from __future__ import annotations

from typing import Any

import numpy as np

from ._validate import as_int_array, recycle
from .rle import Rle


def coverage_k(
    start: np.ndarray,
    width: np.ndarray,
    shift: Any = 0,
    length: int | None = None,
    weight: Any = 1,
) -> Rle:
    """Depth of coverage at positions ``1..length`` as an ``Rle``.

    ``length`` defaults to the largest (shifted) end, floored at 0; ranges are
    clipped to ``[1, length]``. ``weight`` (scalar or per-range) multiplies each
    range's contribution; the result is integer unless a weight is non-integral.
    """
    n = len(start)
    sh = recycle(as_int_array(shift, "shift"), n, "shift") if n else np.array([], dtype=np.int64)
    s = start + sh
    e = s + width - 1
    w = np.asarray(weight)
    if w.ndim == 0:
        w = np.repeat(w, n)
    if len(w) != n:
        raise ValueError("'weight' must be a single value or have the length of the ranges")
    if w.dtype.kind in "iub":
        w = w.astype(np.int64)
    else:
        w = w.astype(np.float64)
        if np.all(w == np.floor(w)):
            w = w.astype(np.int64)
    if length is None:
        length = int(max(0, e.max())) if n else 0
    if length < 0:
        raise ValueError("'width' must be non-negative")
    s = np.maximum(s, 1)
    e = np.minimum(e, length)
    ok = (width > 0) & (e >= s)
    s, e, w = s[ok], e[ok], w[ok]
    if len(s) == 0:
        return Rle(np.zeros(1, dtype=w.dtype), np.array([length])) if length else Rle(
            np.array([], dtype=w.dtype), np.array([], dtype=np.int64)
        )
    pos = np.concatenate((s - 1, e))
    val = np.concatenate((w, -w))
    uniq, inv = np.unique(pos, return_inverse=True)
    delta = np.zeros(len(uniq), dtype=val.dtype)
    np.add.at(delta, inv, val)
    depth = np.cumsum(delta)
    if uniq[0] > 0:
        uniq = np.concatenate(([0], uniq))
        depth = np.concatenate((np.zeros(1, dtype=depth.dtype), depth))
    lengths = np.diff(np.concatenate((uniq, [length])))
    return Rle(depth, lengths)
