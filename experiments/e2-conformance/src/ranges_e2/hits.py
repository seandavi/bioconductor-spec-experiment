"""Hits: pairs (query index, subject index), 0-based, grouped by query."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True)
class Hits:
    query_hits: npt.NDArray[np.int64]
    subject_hits: npt.NDArray[np.int64]
    query_length: int
    subject_length: int

    def __post_init__(self) -> None:
        for name in ("query_hits", "subject_hits"):
            arr = np.array(getattr(self, name), copy=True)
            arr.setflags(write=False)
            object.__setattr__(self, name, arr)

    def __len__(self) -> int:
        return len(self.query_hits)

    def pairs(self) -> list[tuple[int, int]]:
        return list(zip(self.query_hits.tolist(), self.subject_hits.tolist()))
