"""Minimal run-length-encoded vector, used for coverage results."""

from __future__ import annotations

import numpy as np

from ._meta import freeze


class Rle:
    """Run-length encoding: `values[i]` repeated `lengths[i]` times."""

    __slots__ = ("_values", "_lengths")

    def __init__(self, values: np.ndarray, lengths: np.ndarray) -> None:
        v = np.asarray(values)
        n = np.asarray(lengths, dtype=np.int64)
        if v.ndim != 1 or n.ndim != 1 or len(v) != len(n):
            raise ValueError("values and lengths must be 1-d arrays of equal length")
        if (n <= 0).any():
            raise ValueError("lengths must be positive")
        self._values, self._lengths = freeze(v, v.dtype), freeze(n, np.int64)

    @classmethod
    def encode(cls, x: np.ndarray) -> Rle:
        """Run-length encode a plain array (adjacent equal values merge)."""
        x = np.asarray(x)
        if len(x) == 0:
            return cls(x, np.empty(0, dtype=np.int64))
        starts = np.flatnonzero(np.concatenate([[True], x[1:] != x[:-1]]))
        lengths = np.diff(np.append(starts, len(x)))
        return cls(x[starts], lengths)

    @property
    def values(self) -> np.ndarray:
        return self._values

    @property
    def lengths(self) -> np.ndarray:
        return self._lengths

    def __len__(self) -> int:
        return int(self._lengths.sum())

    def decode(self) -> np.ndarray:
        """Expand to a plain array (R: as.vector)."""
        return np.repeat(self._values, self._lengths)

    def __repr__(self) -> str:
        return f"Rle(length={len(self)}, runs={len(self._values)})"
