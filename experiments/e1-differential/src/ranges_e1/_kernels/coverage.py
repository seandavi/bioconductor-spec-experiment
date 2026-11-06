"""Coverage as run-length data, computed from sorted breakpoints (sparse: width may be huge)."""

from __future__ import annotations

import numpy as np


def _plain(
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


def _merge(vals: np.ndarray, lens: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Merge adjacent equal runs (NaN equals NaN)."""
    if len(vals) == 0:
        return vals, lens.astype(np.int64)
    same = (vals[1:] == vals[:-1]) | (np.isnan(vals[1:]) & np.isnan(vals[:-1]))
    idx = np.flatnonzero(np.concatenate([[True], ~same]))
    return vals[idx], np.add.reduceat(lens, idx).astype(np.int64)


def _blank(vals: np.ndarray, lens: np.ndarray, mask_vals: np.ndarray, mask_lens: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Set `vals` to NaN wherever the second run list is positive (both cover the same positions)."""
    cuts = np.union1d(np.cumsum(lens), np.cumsum(mask_lens))
    i = np.searchsorted(np.cumsum(lens), cuts, side="left")
    j = np.searchsorted(np.cumsum(mask_lens), cuts, side="left")
    out = np.where(mask_vals[j] > 0, np.nan, vals[i].astype(np.float64))
    return _merge(out, np.diff(np.concatenate([[0], cuts])))


def coverage(
    start: np.ndarray, width: np.ndarray, weight: np.ndarray, total_width: int
) -> tuple[np.ndarray, np.ndarray]:
    """Weighted coverage of positions 1..total_width as (values, lengths) runs.

    Ranges are clipped to [1, total_width]. `weight` may hold NaN (missing); then the result
    follows R's engine, whose behaviour was reverse-engineered on 4000 random layouts (see
    `_with_missing`).
    """
    if weight.dtype.kind != "f" or not np.isnan(weight).any():
        return _plain(start, width, weight, total_width)
    return _with_missing(start, width, weight, total_width)


def _with_missing(
    start: np.ndarray, width: np.ndarray, weight: np.ndarray, total: int
) -> tuple[np.ndarray, np.ndarray]:
    """R's coverage() with NA weights. Findings, all in input order of the ranges:

    - R picks an internal method: "sort" if total >= 4 * length(x), else "hash".
    - A missing weight only blanks the positions of its own range when the ranges form a
      partition: first start <= 1, each start = previous end + 1, last end >= total.
      Otherwise every position from the first missing range's (clipped) start to `total` is NA.
    - "sort" ignores empty and out-of-window ranges everywhere (including the partition test's
      NA pool but not its contiguity test); "hash" keeps empties, counts missing ranges that lie
      outside the window for the poison start, and tests contiguity on the ranges it kept.
    """
    n = len(start)
    end = start + width - 1
    lo, hi = np.maximum(start, 1), np.minimum(end, total)
    nonempty_in = (width > 0) & (lo <= hi)
    missing = np.isnan(weight)
    sort_method = total >= 4 * n
    if sort_method:
        data = nonempty_in
        order = np.arange(n)  # contiguity is tested on every range, in input order
        pool = data & missing
    else:
        data = (width == 0) | nonempty_in
        order = np.flatnonzero(data)
        pool = missing
    known = data & ~missing
    values, lengths = _plain(start[known], width[known], weight[known], total)
    if not pool.any():
        return values, lengths
    seq_s, seq_e = start[order], end[order]
    partition = (
        len(order) > 0 and seq_s[0] <= 1 and bool((seq_s[1:] == seq_e[:-1] + 1).all()) and seq_e[-1] >= total
    )
    if partition:
        nv, nl = _plain(start[data & missing], width[data & missing], np.ones(int((data & missing).sum())), total)
        return _blank(values, lengths, nv, nl)
    first_na = int(np.maximum(start[pool], 1).min())
    if first_na > total:
        return values, lengths
    return _poison_from(values, lengths, first_na, total)


def _poison_from(vals: np.ndarray, lens: np.ndarray, first_na: int, total: int) -> tuple[np.ndarray, np.ndarray]:
    """Keep positions < first_na, make first_na..total missing."""
    cut = first_na - 1
    cum = np.cumsum(lens)
    keep = int(np.searchsorted(cum, cut, side="left")) + 1 if cut > 0 else 0
    hv, hl = vals[:keep].astype(np.float64), lens[:keep].copy()
    if keep:
        hl[-1] -= cum[keep - 1] - cut
    return _merge(np.concatenate([hv, [np.nan]]), np.concatenate([hl, [total - cut]]))
