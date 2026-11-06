"""Run an IRanges-style kernel per (seqlevel, strand) group of a GRanges.

The group key is ``seqcode * 3 + strand`` (strand 0 +, 1 -, 2 *); with
``ignore_strand`` every range is keyed as ``*``. Groups are visited in key
order, so output is sorted by seqlevel, then strand (+ - *), then position,
as R does. No algorithm is written twice: the kernels in inter.py / setops.py
see only start/width arrays.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

Coords = tuple[np.ndarray, np.ndarray]
# kernel result: (start, width, members) where members index the kernel's input
Result = tuple[np.ndarray, np.ndarray, list[np.ndarray] | None]


def keys(seqcodes: np.ndarray, strand: np.ndarray, ignore_strand: bool) -> np.ndarray:
    return seqcodes * 3 + (2 if ignore_strand else strand.astype(np.int64))


class Collected:
    """Accumulates per-group kernel output into flat arrays."""

    def __init__(self) -> None:
        self.codes: list[np.ndarray] = []
        self.strands: list[np.ndarray] = []
        self.start: list[np.ndarray] = []
        self.width: list[np.ndarray] = []
        self.members: list[np.ndarray] = []
        self.has_members = False

    def add(self, key: int, start: np.ndarray, width: np.ndarray, members: list[np.ndarray] | None = None) -> None:
        n = len(start)
        self.codes.append(np.full(n, key // 3, dtype=np.int64))
        self.strands.append(np.full(n, key % 3, dtype=np.int8))
        self.start.append(start)
        self.width.append(width)
        if members is not None:
            self.has_members = True
            self.members.extend(members)

    def flat(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        def cat(parts: list[np.ndarray], dtype: type) -> np.ndarray:
            return np.concatenate(parts) if parts else np.zeros(0, dtype=dtype)

        return cat(self.codes, np.int64), cat(self.strands, np.int8), cat(self.start, np.int64), cat(self.width, np.int64)


def by_group(
    seqcodes: np.ndarray,
    strand: np.ndarray,
    start: np.ndarray,
    width: np.ndarray,
    kernel: Callable[[np.ndarray, np.ndarray], Result],
    ignore_strand: bool,
) -> Collected:
    """Apply ``kernel`` to each non-empty group."""
    k = keys(seqcodes, strand, ignore_strand)
    out = Collected()
    for key in sorted(set(k.tolist())):
        idx = np.flatnonzero(k == key)
        s, w, members = kernel(start[idx], width[idx])
        out.add(key, s, w, None if members is None else [idx[m] for m in members])
    return out


def by_group_pair(
    xk: np.ndarray,
    x: Coords,
    yk: np.ndarray,
    y: Coords,
    kernel: Callable[[Coords, Coords], Coords],
) -> Collected:
    out = Collected()
    for key in sorted(set(xk.tolist()) | set(yk.tolist())):
        xi, yi = np.flatnonzero(xk == key), np.flatnonzero(yk == key)
        s, w = kernel((x[0][xi], x[1][xi]), (y[0][yi], y[1][yi]))
        out.add(key, s, w)
    return out
