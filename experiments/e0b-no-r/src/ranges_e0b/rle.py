"""Minimal run-length encoding, enough to return coverage."""

from typing import Any

import numpy as np
import numpy.typing as npt

from ._sew import IntArray


class Rle:
    """Run-length encoded vector: `values[i]` repeated `lengths[i]` times."""

    __slots__ = ("_values", "_lengths")

    def __init__(self, values: Any, lengths: Any) -> None:
        v = np.asarray(values)
        n = np.asarray(lengths, dtype=np.int64)
        if v.ndim != 1 or v.shape != n.shape:
            raise ValueError("values and lengths must be 1-d and the same length")
        if np.any(n <= 0):
            raise ValueError("run lengths must be positive")
        self._values, self._lengths = v, n

    @classmethod
    def from_array(cls, x: Any) -> "Rle":
        """Encode a plain array (adjacent equal values are merged)."""
        a = np.asarray(x)
        if len(a) == 0:
            return cls(a, np.empty(0, np.int64))
        starts = np.flatnonzero(np.concatenate(([True], a[1:] != a[:-1])))
        return cls(a[starts], np.diff(np.append(starts, len(a))))

    @property
    def values(self) -> npt.NDArray[Any]:
        return self._values

    @property
    def lengths(self) -> IntArray:
        return self._lengths

    runvalue = values
    runlength = lengths

    def __len__(self) -> int:
        return int(self._lengths.sum())

    def to_numpy(self) -> npt.NDArray[Any]:
        return np.repeat(self._values, self._lengths)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Rle):
            return NotImplemented
        a, b = Rle.from_array(self.to_numpy()), Rle.from_array(other.to_numpy())
        return bool(np.array_equal(a._values, b._values) and np.array_equal(a._lengths, b._lengths))

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        runs = ", ".join(f"{v}x{n}" for v, n in zip(self._values[:8].tolist(), self._lengths[:8].tolist()))
        more = ", ..." if len(self._values) > 8 else ""
        return f"Rle(length={len(self)}, runs={len(self._values)}: {runs}{more})"
