"""``Rle``: minimal run-length vector, the result type of ``coverage`` (COV.3)."""

from __future__ import annotations

import numpy as np

from ranges_e3._types import IntArray


class Rle:
    """Runs of equal int64 values; ``lengths > 0``, adjacent runs differ."""

    __slots__ = ("_lengths", "_values")

    def __init__(self, lengths: IntArray, values: IntArray) -> None:
        if len(lengths) != len(values):
            raise ValueError("lengths and values must have equal length")
        if (np.asarray(lengths) <= 0).any():
            raise ValueError("run lengths must be positive")
        self._lengths = np.asarray(lengths, dtype=np.int64)
        self._values = np.asarray(values, dtype=np.int64)

    @property
    def lengths(self) -> IntArray:
        return self._lengths.copy()

    @property
    def values(self) -> IntArray:
        return self._values.copy()

    def __len__(self) -> int:
        return int(self._lengths.sum())

    def to_numpy(self) -> IntArray:
        """Decoded values, one per position."""
        return np.repeat(self._values, self._lengths)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Rle):
            return NotImplemented
        return bool(
            np.array_equal(self._lengths, other._lengths) and np.array_equal(self._values, other._values)
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        runs = ", ".join(f"{n}x{v}" for n, v in zip(self._lengths, self._values))
        return f"Rle({runs})"
