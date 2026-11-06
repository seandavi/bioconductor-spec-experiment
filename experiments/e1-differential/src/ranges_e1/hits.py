"""Hits: sparse query/subject matches, sorted by query (0-based indices)."""

from __future__ import annotations

import numpy as np

from ._meta import freeze


class Hits:
    """Query/subject index pairs. Indices are 0-based (ADR-0001)."""

    __slots__ = ("_q", "_s", "_nq", "_ns")

    def __init__(self, query_hits: np.ndarray, subject_hits: np.ndarray, query_length: int, subject_length: int) -> None:
        q = np.asarray(query_hits, dtype=np.int64)
        s = np.asarray(subject_hits, dtype=np.int64)
        if len(q) != len(s):
            raise ValueError("query_hits and subject_hits must have the same length")
        if len(q) and (q.min() < 0 or q.max() >= query_length or s.min() < 0 or s.max() >= subject_length):
            raise ValueError("hit index out of range")
        order = np.argsort(q, kind="stable")  # R sorts by query only; subject order is the engine's
        q, s = q[order], s[order]
        q, s = freeze(q, np.int64), freeze(s, np.int64)
        self._q, self._s, self._nq, self._ns = q, s, int(query_length), int(subject_length)

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

    def __repr__(self) -> str:
        return f"Hits(n={len(self)}, query_length={self._nq}, subject_length={self._ns})"

    def count_by_query(self) -> np.ndarray:
        """Number of hits per query (R: countQueryHits)."""
        return np.bincount(self._q, minlength=self._nq).astype(np.int64)

    def select(self, how: str) -> np.ma.MaskedArray:
        """Pick one subject index per query: "first" (min), "last" (max), or "arbitrary" (R's: the query's first hit in engine order)."""
        if how not in ("first", "last", "arbitrary"):
            raise ValueError("'arg' should be one of \"all\", \"first\", \"last\", \"arbitrary\"")
        has = np.zeros(self._nq, dtype=bool)
        has[self._q] = True
        if how == "arbitrary":
            out = np.zeros(self._nq, dtype=np.int64)
            out[self._q[::-1]] = self._s[::-1]  # first write in stored order wins
        elif how == "last":
            out = np.full(self._nq, -1, dtype=np.int64)
            np.maximum.at(out, self._q, self._s)
        else:
            out = np.full(self._nq, np.iinfo(np.int64).max, dtype=np.int64)
            np.minimum.at(out, self._q, self._s)
        return np.ma.masked_array(np.where(has, out, 0), mask=~has)


def overlap_result(
    q: np.ndarray,
    s: np.ndarray,
    n_query: int,
    n_subject: int,
    *,
    select: str,
    self_hits: bool,
    drop_self: bool = False,
    drop_redundant: bool = False,
    query_is_shorter: bool | np.ndarray = False,
) -> Hits | np.ma.MaskedArray:
    """Turn raw overlap pairs into R's findOverlaps result.

    `select="arbitrary"` mimics R's engines: with self-overlaps or a shorter query R reports the
    last-seen subject (max index); otherwise the first hit in engine order.
    """
    if select not in ("all", "first", "last", "arbitrary"):
        raise ValueError("'arg' should be one of \"all\", \"first\", \"last\", \"arbitrary\"")
    if self_hits:
        keep = np.ones(len(q), dtype=bool)
        if drop_self:
            keep &= q != s
        if drop_redundant:  # drop (i, j), i > j, when its mirror (j, i) is also a hit
            key = set((q * n_query + s).tolist())
            mirror = np.array([(b * n_query + a) in key for a, b in zip(q.tolist(), s.tolist())], dtype=bool)
            keep &= ~((q > s) & mirror)
        q, s = q[keep], s[keep]
    hits = Hits(q, s, n_query, n_subject)
    if select == "all":
        return hits
    if select == "arbitrary":
        if self_hits:
            return hits.select("last")
        shorter = np.broadcast_to(np.asarray(query_is_shorter, dtype=bool), (n_query,))
        if shorter.any():  # per-query engine mode (GRanges runs it per seqlevel)
            return np.ma.where(shorter, hits.select("last"), hits.select("arbitrary"))
    return hits.select(select)
