"""Coverage as run-length data, computed from sorted breakpoints (sparse: width may be huge)."""

from __future__ import annotations

import numpy as np


def coverage(
    start: np.ndarray, width: np.ndarray, weight: np.ndarray, total_width: int
) -> tuple[np.ndarray, np.ndarray]:
    """Sum of weights covering positions 1..total_width; returns (values, lengths) runs.

    Ranges are clipped to [1, total_width]; empty ranges contribute nothing.
    """
    end = start + width - 1
    a = np.maximum(start, 1)
    b = np.minimum(end, total_width)
    ok = (width > 0) & (a <= b)
    a, b, w = a[ok], b[ok], weight[ok]
    pos = np.concatenate([a, b + 1])
    delta = np.concatenate([w, -w])
    upos, inv = np.unique(pos, return_inverse=True)
    step = np.zeros(len(upos), dtype=delta.dtype)
    np.add.at(step, inv, delta)
    level = np.cumsum(step)
    # segments: [1, upos[0]) at 0, then [upos[i], upos[i+1]) at level[i]; last ends at total_width+1
    bounds = np.concatenate([[1], upos, [total_width + 1]])
    vals = np.concatenate([np.zeros(1, dtype=delta.dtype), level])
    lens = np.diff(bounds)
    keep = lens > 0
    vals, lens = vals[keep], lens[keep]
    if len(vals) == 0:
        return vals, lens.astype(np.int64)
    first = np.concatenate([[True], vals[1:] != vals[:-1]])
    idx = np.flatnonzero(first)
    return vals[idx], np.add.reduceat(lens, idx).astype(np.int64)
