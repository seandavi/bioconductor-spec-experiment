"""``Hits``: query/subject index pairs (OVL.9)."""

from __future__ import annotations

import numpy as np

from ranges_e3._types import IntArray


class Hits:
    """Hit pairs ordered by (query, subject); indices are 0-based (CONV.3)."""

    __slots__ = ("_query", "_subject", "_nq", "_ns")

    def __init__(self, query_hits: IntArray, subject_hits: IntArray, query_length: int, subject_length: int) -> None:
        if len(query_hits) != len(subject_hits):
            raise ValueError("query_hits and subject_hits must have equal length")
        self._query = np.asarray(query_hits, dtype=np.int64)
        self._subject = np.asarray(subject_hits, dtype=np.int64)
        self._nq, self._ns = int(query_length), int(subject_length)

    @property
    def query_hits(self) -> IntArray:
        return self._query.copy()

    @property
    def subject_hits(self) -> IntArray:
        return self._subject.copy()

    @property
    def query_length(self) -> int:
        return self._nq

    @property
    def subject_length(self) -> int:
        return self._ns

    def __len__(self) -> int:
        return len(self._query)

    def pairs(self) -> list[tuple[int, int]]:
        return list(zip(self._query.tolist(), self._subject.tolist()))

    def __repr__(self) -> str:
        return f"Hits(n={len(self)}, query_length={self._nq}, subject_length={self._ns})"
