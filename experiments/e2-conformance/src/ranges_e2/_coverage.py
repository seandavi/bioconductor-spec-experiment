"""Coverage kernel: sparse event sweep producing run-length output."""

from __future__ import annotations

import numpy as np

from ._util import IntArray, rep_len
from .rle import Rle


def coverage_rle(
    start: IntArray,
    width: IntArray,
    length: int,
    shift: IntArray,
    weight: np.ndarray,
) -> Rle:
    """Weighted coverage over positions 1..length. Positions < 1 and > length are dropped.

    weight int dtype -> int32 values; float dtype -> float64 values.
    """
    is_int = weight.dtype.kind in "iu"
    dtype = np.int64 if is_int else np.float64
    out_dtype = np.int32 if is_int else np.float64
    if length < 0:
        raise ValueError("width cannot be negative")
    if length == 0:
        return Rle(np.zeros(0, dtype=out_dtype), np.zeros(0, dtype=np.int64))
    n = len(start)
    s = start + rep_len(shift, n, "shift") if n else start
    w = rep_len(weight, n, "weight").astype(dtype) if n else weight.astype(dtype)
    lo = np.maximum(s, 1)
    hi = np.minimum(s + width, length + 1)
    keep = (width > 0) & (hi > lo)
    lo, hi, w = lo[keep], hi[keep], w[keep]
    if len(lo) == 0:
        return Rle(np.zeros(1, dtype=out_dtype), np.array([length], dtype=np.int64))
    pos = np.concatenate((lo, hi))
    delta = np.concatenate((w, -w))
    uniq, inv = np.unique(pos, return_inverse=True)
    agg = np.zeros(len(uniq), dtype=dtype)
    np.add.at(agg, inv, delta)
    level = np.cumsum(agg)  # value on [uniq[k], uniq[k+1])
    bounds = np.concatenate(([1], uniq, [length + 1]))
    values = np.concatenate(([0], level)).astype(dtype)  # value on [bounds[k], bounds[k+1])
    lengths = np.diff(bounds)
    keep_run = lengths > 0
    values, lengths = values[keep_run], lengths[keep_run]
    change = np.ones(len(values), dtype=bool)
    change[1:] = values[1:] != values[:-1]
    group = np.cumsum(change) - 1
    merged_len = np.zeros(group[-1] + 1, dtype=np.int64)
    np.add.at(merged_len, group, lengths)
    return Rle(values[change].astype(out_dtype), merged_len)
