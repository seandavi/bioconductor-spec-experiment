"""Coverage kernel: per-position depth as a dense array, via a difference array."""

from __future__ import annotations

from typing import Any

import numpy as np

from ._arrays import as_int
from .intra import arg_vector


def coverage_dense(
    start: np.ndarray, width: np.ndarray, shift: Any, length: int | None, weight: Any
) -> np.ndarray:
    """Depth at positions 1..L. L is ``length`` or max(end + shift) (zero-width ranges count)."""
    n = len(start)
    sh = arg_vector(shift, n, "shift") if n else np.zeros(0, np.int64)
    w = np.asarray(weight)
    if w.ndim == 0 or len(w) == 1:
        w = np.repeat(w.reshape(-1)[:1], n)
    elif len(w) != n:
        raise ValueError(f"'weight' has length {len(w)}, expected 1 or {n}")
    if w.dtype.kind not in "iuf":
        raise TypeError("'weight' must be numeric")
    dtype = np.float64 if w.dtype.kind == "f" else np.int64
    s = start + sh
    e = s + width - 1
    total = int(max(e.max(), 0)) if n else 0
    if length is not None:
        if length < 0:
            raise ValueError("'width' must be non-negative")
        total = length
    delta = np.zeros(total + 2, dtype=dtype)
    keep = np.flatnonzero(width > 0)
    lo = np.maximum(s[keep], 1)
    hi = np.minimum(e[keep], total)
    ok = hi >= lo
    np.add.at(delta, lo[ok], w[keep][ok].astype(dtype))
    np.add.at(delta, hi[ok] + 1, -w[keep][ok].astype(dtype))
    return np.cumsum(delta)[1 : total + 1]


def check_length(width: Any) -> int | None:
    if width is None:
        return None
    return int(as_int(width, "width")[0])
