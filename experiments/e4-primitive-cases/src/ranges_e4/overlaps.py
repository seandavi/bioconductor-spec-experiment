"""Overlap kernel over start/width arrays; everything else builds on ``overlap_pairs``.

Predicates (q = query, s = subject, ends inclusive, g = maxgap):
  any     qs - se - 1 <= g  and  ss - qe - 1 <= g
  start   |qs - ss| <= g          end    |qe - se| <= g
  equal   both of the above       within ss <= qs, qe <= se, and
                                          (g > 0 -> (se-ss) - (qe-qs) <= g)
``maxgap=-1`` means 0 for start/end/equal and "no limit" for within.
``minoverlap > 0`` additionally requires ``min(qe,se) - max(qs,ss) + 1 >= minoverlap``.
Probed against R, including zero-width ranges and the within/maxgap=0 quirk.
Candidates are bounded by a sorted-start prefix per query; the filter is exact.
"""

from __future__ import annotations

import numpy as np

from .hits import Hits

TYPES = ("any", "start", "end", "within", "equal")
SELECTS = ("all", "first", "last", "arbitrary")


def check_args(maxgap: int, minoverlap: int, type_: str) -> None:
    if type_ not in TYPES:
        raise ValueError(f"'type' must be one of {TYPES}")
    if maxgap < -1:
        raise ValueError("'maxgap' must be >= -1")
    if minoverlap < 0:
        raise ValueError("'minoverlap' cannot be negative")
    if type_ == "any" and maxgap != -1 and minoverlap != 0:
        raise ValueError("for type 'any', at least one of 'maxgap' and 'minoverlap' must be default")


def overlap_pairs(
    qs: np.ndarray, qw: np.ndarray, ss: np.ndarray, sw: np.ndarray, maxgap: int, minoverlap: int, type_: str
) -> tuple[np.ndarray, np.ndarray]:
    """All (query, subject) positions satisfying the predicate, sorted by query then subject."""
    check_args(maxgap, minoverlap, type_)
    qe, se = qs + qw - 1, ss + sw - 1
    order = np.argsort(ss, kind="stable")
    ss_sorted = ss[order]
    reach = max(maxgap, 0) + 1
    hi = np.searchsorted(ss_sorted, qe + reach, side="right")
    q_out: list[np.ndarray] = []
    s_out: list[np.ndarray] = []
    for i in range(len(qs)):
        cand = order[: hi[i]]
        if not len(cand):
            continue
        cs, ce = ss[cand], se[cand]
        a, b = qs[i], qe[i]
        match type_:
            case "any":
                ok = (a - ce - 1 <= maxgap) & (cs - b - 1 <= maxgap)
            case "start":
                ok = np.abs(a - cs) <= max(maxgap, 0)
            case "end":
                ok = np.abs(b - ce) <= max(maxgap, 0)
            case "equal":
                g = max(maxgap, 0)
                ok = (np.abs(a - cs) <= g) & (np.abs(b - ce) <= g)
            case _:  # within
                ok = (cs <= a) & (b <= ce)
                if maxgap > 0:
                    ok &= (ce - cs) - (b - a) <= maxgap
        if minoverlap > 0:
            ok &= np.minimum(b, ce) - np.maximum(a, cs) + 1 >= minoverlap
        hit = np.sort(cand[ok])
        if len(hit):
            q_out.append(np.full(len(hit), i, dtype=np.int64))
            s_out.append(hit.astype(np.int64))
    if not q_out:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(q_out), np.concatenate(s_out)


def select_hits(q: np.ndarray, s: np.ndarray, n_query: int, select: str) -> np.ndarray:
    """Collapse sorted pairs to one subject index per query (masked where no hit)."""
    if select not in SELECTS:
        raise ValueError(f"'select' must be one of {SELECTS}")
    out = np.zeros(n_query, dtype=np.int64)
    has = np.zeros(n_query, dtype=bool)
    if len(q):
        # pairs are sorted by (query, subject): first = first row per query, last = last row
        pick_last = select == "last"
        idx = np.flatnonzero(np.r_[q[1:] != q[:-1], True]) if pick_last else np.flatnonzero(np.r_[True, q[1:] != q[:-1]])
        out[q[idx]] = s[idx]
        has[q[idx]] = True
    return np.ma.MaskedArray(out, mask=~has)


def to_hits(q: np.ndarray, s: np.ndarray, n_query: int, n_subject: int) -> Hits:
    return Hits(q, s, n_query, n_subject)
