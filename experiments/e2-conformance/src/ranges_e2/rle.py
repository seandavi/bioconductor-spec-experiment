"""Minimal Factor and Rle (run-length vector) used for coverage and GRanges accessors."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True)
class Factor:
    """Categorical vector: 0-based `codes` into `levels`; -1 means NA."""

    levels: tuple[str, ...]
    codes: npt.NDArray[np.int64]

    def __post_init__(self) -> None:
        codes = np.array(self.codes, copy=True)
        codes.setflags(write=False)
        object.__setattr__(self, "codes", codes)

    def __len__(self) -> int:
        return len(self.codes)

    def to_list(self) -> list[str | None]:
        return [None if c < 0 else self.levels[c] for c in self.codes.tolist()]


@dataclass(frozen=True)
class Rle:
    """Run-length encoded vector. `values` is a numpy array or a Factor (one entry per run)."""

    values: npt.NDArray[np.generic] | Factor
    lengths: npt.NDArray[np.int64]

    def __post_init__(self) -> None:
        lengths = np.array(self.lengths, copy=True)
        lengths.setflags(write=False)
        object.__setattr__(self, "lengths", lengths)
        if isinstance(self.values, np.ndarray):
            values = np.array(self.values, copy=True)
            values.setflags(write=False)
            object.__setattr__(self, "values", values)

    def __len__(self) -> int:
        return int(self.lengths.sum())

    @classmethod
    def from_dense(cls, dense: npt.NDArray[np.generic]) -> Rle:
        """Compress a dense array into runs of equal values."""
        if len(dense) == 0:
            return cls(dense[:0], np.zeros(0, dtype=np.int64))
        change = np.nonzero(dense[1:] != dense[:-1])[0] + 1
        starts = np.concatenate(([0], change))
        lengths = np.diff(np.concatenate((starts, [len(dense)]))).astype(np.int64)
        return cls(dense[starts], lengths)

    @classmethod
    def from_factor(cls, factor: Factor) -> Rle:
        rle = cls.from_dense(factor.codes)
        assert isinstance(rle.values, np.ndarray)
        return cls(Factor(factor.levels, rle.values.astype(np.int64)), rle.lengths)

    def to_dense(self) -> npt.NDArray[np.generic]:
        if isinstance(self.values, Factor):
            return np.repeat(self.values.codes, self.lengths)
        return np.repeat(self.values, self.lengths)
