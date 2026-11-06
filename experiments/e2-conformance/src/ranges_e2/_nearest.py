"""nearest / precede / follow / distance kernels. Per-query vectorised over subjects.

Boundaries are plain arrays so a GRanges layer can pass strand/seqname masks per query.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import numpy.typing as npt

from ._util import IntArray, rep_len

BoolArray = npt.NDArray[np.bool_]
# mask(query_index) -> bool array over subjects: True where the pair is eligible
MaskFn = Callable[[int], BoolArray]
_INF = np.iinfo(np.int64).max


def pair_distance(qs: IntArray, qe: IntArray, ss: IntArray, se: IntArray) -> IntArray:
    """Number of positions between two ranges; 0 if they overlap or are adjacent."""
    return np.maximum(0, np.maximum(qs, ss) - np.minimum(qe, se) - 1)


def nearest_ties(
    qs: IntArray, qw: IntArray, ss: IntArray, sw: IntArray, mask: MaskFn | None = None, self_exclude: bool = False
) -> list[IntArray]:
    """Per query: all eligible subjects at the minimum distance (empty array if none)."""
    qe = qs + qw - 1
    se = ss + sw - 1
    out: list[IntArray] = []
    for i in range(len(qs)):
        dist = pair_distance(qs[i], qe[i], ss, se).astype(np.int64)
        if mask is not None:
            dist = np.where(mask(i), dist, _INF)
        if self_exclude:
            dist[i] = _INF
        if len(dist) == 0 or dist.min() == _INF:
            out.append(np.zeros(0, dtype=np.int64))
        else:
            out.append(np.nonzero(dist == dist.min())[0].astype(np.int64))
    return out


def directional_ties(
    qs: IntArray,
    qw: IntArray,
    ss: IntArray,
    sw: IntArray,
    kind: str,
    reverse: MaskFn | None = None,
    mask: MaskFn | None = None,
    self_exclude: bool = False,
) -> list[IntArray]:
    """precede/follow candidates per query: eligible subjects at the minimum directional distance.

    precede (forward): subject starts after end(q), distance = start(s) - end(q).
    follow  (forward): subject ends before start(q), distance = start(q) - end(s).
    `reverse(i)` marks subjects for which the direction is flipped (GRanges '-' strand).
    """
    if kind not in ("precede", "follow"):
        raise ValueError(kind)
    qe = qs + qw - 1
    se = ss + sw - 1
    out: list[IntArray] = []
    for i in range(len(qs)):
        after = np.where(ss > qe[i], ss - qe[i], _INF)  # subject lies after q
        before = np.where(se < qs[i], qs[i] - se, _INF)  # subject lies before q
        fwd, rev = (after, before) if kind == "precede" else (before, after)
        if reverse is None:
            dist = fwd
        else:
            dist = np.where(reverse(i), rev, fwd)
        if mask is not None:
            dist = np.where(mask(i), dist, _INF)
        if self_exclude:
            dist = dist.copy()
            dist[i] = _INF
        if len(dist) == 0 or dist.min() == _INF:
            out.append(np.zeros(0, dtype=np.int64))
        else:
            out.append(np.nonzero(dist == dist.min())[0].astype(np.int64))
    return out


def distance_ranges(xs: IntArray, xw: IntArray, ys: IntArray, yw: IntArray) -> IntArray:
    """Parallel distance; the shorter argument recycles to the longer length."""
    n = max(len(xs), len(ys)) if len(xs) and len(ys) else 0
    xs, xw, ys, yw = (rep_len(a, n) for a in (xs, xw, ys, yw))
    return pair_distance(xs, xs + xw - 1, ys, ys + yw - 1)
