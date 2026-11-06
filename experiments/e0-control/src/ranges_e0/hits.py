"""Hits: the result of ``find_overlaps(select="all")``."""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np


class Hits:
    """Query/subject index pairs (0-based), sorted by query then subject."""

    __slots__ = ("_q", "_s", "_nq", "_ns")

    def __init__(
        self,
        query_hits: np.ndarray,
        subject_hits: np.ndarray,
        query_length: int,
        subject_length: int,
    ) -> None:
        q = np.asarray(query_hits, dtype=np.int64)
        s = np.asarray(subject_hits, dtype=np.int64)
        if q.shape != s.shape or q.ndim != 1:
            raise ValueError("query_hits and subject_hits must be 1-d and of equal length")
        if len(q) and (q.min() < 0 or q.max() >= query_length or s.min() < 0 or s.max() >= subject_length):
            raise ValueError("hit indices out of range")
        order = np.lexsort((s, q))
        self._q = q[order]
        self._s = s[order]
        self._nq = int(query_length)
        self._ns = int(subject_length)

    @property
    def query_hits(self) -> np.ndarray:
        return self._q

    @property
    def subject_hits(self) -> np.ndarray:
        return self._s

    @property
    def query_length(self) -> int:
        return self._nq

    @property
    def subject_length(self) -> int:
        return self._ns

    def __len__(self) -> int:
        return len(self._q)

    def __iter__(self) -> Iterator[tuple[int, int]]:
        return iter(zip(self._q.tolist(), self._s.tolist()))

    def to_numpy(self) -> np.ndarray:
        """An ``(n, 2)`` array of ``(query, subject)`` pairs."""
        return np.column_stack((self._q, self._s))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Hits):
            return NotImplemented
        return (
            self._nq == other._nq
            and self._ns == other._ns
            and np.array_equal(self._q, other._q)
            and np.array_equal(self._s, other._s)
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return (
            f"Hits(n={len(self)}, query_length={self._nq}, subject_length={self._ns}, "
            f"pairs={self.to_numpy().tolist()})"
        )
