"""Run-length vector; only used as the result of ``coverage``."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Rle:
    values: np.ndarray
    lengths: np.ndarray

    def __len__(self) -> int:
        return int(self.lengths.sum())

    @classmethod
    def encode(cls, x: np.ndarray) -> Rle:
        """Run-length encode a dense vector (adjacent equal values merge)."""
        if not len(x):
            return cls(x[:0].copy(), np.zeros(0, np.int64))
        starts = np.flatnonzero(np.r_[True, x[1:] != x[:-1]])
        lengths = np.diff(np.r_[starts, len(x)]).astype(np.int64)
        return cls(x[starts], lengths)

    def decode(self) -> np.ndarray:
        return np.repeat(self.values, self.lengths)
