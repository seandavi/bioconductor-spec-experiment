"""Hits: (query, subject) index pairs with the two vector lengths. 0-based (ADR-0001)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Hits:
    """Pairs sorted by query then subject."""

    query: np.ndarray
    subject: np.ndarray
    query_length: int
    subject_length: int

    def __len__(self) -> int:
        return len(self.query)
