"""Minimal run-length vector, used as the return type of ``coverage``."""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import numpy as np


class Rle:
    """Run-length encoded 1-d vector (``lengths`` and ``values`` of equal size).

    Deliberately tiny: no arithmetic, no slicing. Adjacent runs with equal
    values are always merged (the canonical form R's ``Rle`` also keeps).
    """

    __slots__ = ("_lengths", "_values")

    def __init__(self, values: Any, lengths: Any) -> None:
        values = np.asarray(values)
        lengths = np.asarray(lengths, dtype=np.int64)
        if values.ndim != 1 or lengths.ndim != 1 or len(values) != len(lengths):
            raise ValueError("'values' and 'lengths' must be 1-d and of equal length")
        if np.any(lengths <= 0):
            keep = lengths > 0
            if np.any(lengths < 0):
                raise ValueError("'lengths' must be non-negative")
            values, lengths = values[keep], lengths[keep]
        if len(values) > 1:
            starts = np.concatenate(([True], values[1:] != values[:-1]))
            if not starts.all():
                grp = np.cumsum(starts) - 1
                lengths = np.bincount(grp, weights=lengths).astype(np.int64)
                values = values[starts]
        self._values = values
        self._lengths = lengths

    @classmethod
    def from_array(cls, x: Any) -> Rle:
        """Run-length encode a plain 1-d array."""
        x = np.asarray(x)
        if x.size == 0:
            return cls(x, np.array([], dtype=np.int64))
        starts = np.flatnonzero(np.concatenate(([True], x[1:] != x[:-1])))
        lengths = np.diff(np.concatenate((starts, [len(x)])))
        return cls(x[starts], lengths)

    @property
    def values(self) -> np.ndarray:
        return self._values

    @property
    def lengths(self) -> np.ndarray:
        return self._lengths

    def __len__(self) -> int:
        return int(self._lengths.sum())

    def to_numpy(self) -> np.ndarray:
        """Decompress to a plain array."""
        return np.repeat(self._values, self._lengths)

    def runs(self) -> Iterator[tuple[Any, int]]:
        """Iterate ``(value, run_length)`` pairs."""
        for v, n in zip(self._values.tolist(), self._lengths.tolist()):
            yield v, n

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Rle):
            return NotImplemented
        return bool(
            np.array_equal(self._lengths, other._lengths)
            and np.array_equal(self._values, other._values)
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        body = ", ".join(f"{v}x{n}" for v, n in self.runs())
        return f"Rle(length={len(self)}, runs=[{body}])"
