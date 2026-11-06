"""Hits and the overlap kernel."""

from dataclasses import dataclass

import numpy as np

from ._sew import NO_HIT, BoolArray, IntArray, empty_int

TYPES = ("any", "start", "end", "within", "equal")
_CHUNK = 4_000_000  # max candidate pairs materialised at once


@dataclass(frozen=True, eq=False)
class Hits:
    """Overlap hits: 0-based (query, subject) pairs sorted by query, then subject."""

    query_hits: IntArray
    subject_hits: IntArray
    query_length: int
    subject_length: int

    def __len__(self) -> int:
        return len(self.query_hits)

    def pairs(self) -> list[tuple[int, int]]:
        return list(zip(self.query_hits.tolist(), self.subject_hits.tolist()))

    def equals(self, other: "Hits") -> bool:
        return (
            self.query_length == other.query_length
            and self.subject_length == other.subject_length
            and bool(np.array_equal(self.query_hits, other.query_hits))
            and bool(np.array_equal(self.subject_hits, other.subject_hits))
        )

    def count_per_query(self) -> IntArray:
        return np.bincount(self.query_hits, minlength=self.query_length).astype(np.int64)

    def __repr__(self) -> str:
        return f"Hits({len(self)} hits; {self.query_length} queries x {self.subject_length} subjects)"


def check_args(type: str, maxgap: int, minoverlap: int) -> None:
    if type not in TYPES:
        raise ValueError(f"type must be one of {TYPES}, got {type!r}")
    if maxgap < -1:
        raise ValueError("maxgap must be >= -1")
    if minoverlap < 0:
        raise ValueError("minoverlap must be >= 0")
    if type == "any" and maxgap != -1 and minoverlap != 0:
        raise ValueError("for type='any', maxgap and minoverlap cannot both be non-default")


def hit_mask(
    qs: IntArray, qe: IntArray, ss: IntArray, se: IntArray, type: str, maxgap: int, minoverlap: int
) -> BoolArray:
    """Elementwise overlap predicate for aligned (query, subject) pairs."""
    m = max(maxgap, 0)
    if type == "any":
        q_empty, s_empty = qe < qs, se < ss
        inside = (q_empty & ~s_empty & (ss < qs) & (qs <= se)) | (s_empty & ~q_empty & (qs < ss) & (ss <= qe))
        gap = np.where(inside, -1, np.maximum(qs, ss) - np.minimum(qe, se) - 1)
        ok = gap <= maxgap
    elif type == "start":
        ok = np.abs(qs - ss) <= m
    elif type == "end":
        ok = np.abs(qe - se) <= m
    elif type == "equal":
        ok = (np.abs(qs - ss) <= m) & (np.abs(qe - se) <= m)
    else:  # within
        ok = (qs >= ss) & (qe <= se)
        if maxgap >= 0:
            ok &= (qs - ss) + (se - qe) <= maxgap
    if minoverlap > 0:
        ok &= (np.minimum(qe, se) - np.maximum(qs, ss) + 1) >= minoverlap
    return ok


def _window(
    qs: IntArray, qe: IntArray, type: str, maxgap: int, max_sw: int
) -> tuple[IntArray, IntArray]:
    """Inclusive bounds on subject start that contain every possible hit."""
    m = max(maxgap, 0)
    if type == "any":
        return qs - m - max_sw - 1, qe + m + 1
    if type in ("start", "equal"):
        return qs - m, qs + m
    if type == "end":
        return qe - m - max_sw - 1, qe + m + 1
    return qe - max_sw - 1, qs  # within


def overlap_pairs(
    qs: IntArray,
    qe: IntArray,
    ss: IntArray,
    se: IntArray,
    type: str = "any",
    maxgap: int = -1,
    minoverlap: int = 0,
) -> tuple[IntArray, IntArray]:
    """All hit pairs (query index, subject index), sorted by (query, subject)."""
    check_args(type, maxgap, minoverlap)
    nq, ns = len(qs), len(ss)
    if nq == 0 or ns == 0:
        return empty_int(), empty_int()
    order = np.argsort(ss, kind="stable")
    s_sorted = ss[order]
    max_sw = int(max((se - ss + 1).max(), 0))
    lo_v, hi_v = _window(qs, qe, type, maxgap, max_sw)
    lo = np.searchsorted(s_sorted, lo_v, side="left")
    hi = np.searchsorted(s_sorted, hi_v, side="right")
    counts = np.maximum(hi - lo, 0)
    cum = np.cumsum(counts)
    out_q: list[IntArray] = []
    out_s: list[IntArray] = []
    start_q = 0
    while start_q < nq:
        base = cum[start_q - 1] if start_q else 0
        stop_q = int(np.searchsorted(cum, base + _CHUNK, side="right"))
        stop_q = min(max(stop_q, start_q + 1), nq)
        c = counts[start_q:stop_q]
        total = int(c.sum())
        if total:
            q_idx = np.repeat(np.arange(start_q, stop_q, dtype=np.int64), c)
            offs = np.arange(total, dtype=np.int64) - np.repeat(np.cumsum(c) - c, c)
            s_idx = order[np.repeat(lo[start_q:stop_q], c) + offs]
            ok = hit_mask(qs[q_idx], qe[q_idx], ss[s_idx], se[s_idx], type, maxgap, minoverlap)
            out_q.append(q_idx[ok])
            out_s.append(s_idx[ok])
        start_q = stop_q
    if not out_q:
        return empty_int(), empty_int()
    qi, si = np.concatenate(out_q), np.concatenate(out_s)
    o = np.lexsort((si, qi))
    return qi[o], si[o]


def select_hits(qi: IntArray, si: IntArray, n_query: int, select: str) -> IntArray:
    """Collapse sorted pairs to one subject per query (NO_HIT if none)."""
    if select not in ("first", "last", "arbitrary"):
        raise ValueError(f"select must be 'all', 'first', 'last' or 'arbitrary', got {select!r}")
    out = np.full(n_query, NO_HIT, dtype=np.int64)
    if select == "last":
        out[qi] = si
    else:  # arbitrary is implemented as first; deterministic
        out[qi[::-1]] = si[::-1]
    return out
