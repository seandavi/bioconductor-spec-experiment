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

    def __len__(self) -> int:
        return len(self.query_hits)

    def pairs(self) -> list[tuple[int, int]]:
        return list(zip(self.query_hits.tolist(), self.subject_hits.tolist()))
