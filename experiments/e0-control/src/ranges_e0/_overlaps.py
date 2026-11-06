"""Overlap engine.

``overlap_pairs`` is the single kernel behind find/count/any/subset overlaps.
Strategy: sort the subject by start; for each query a ``searchsorted`` gives a
contiguous block of candidate subjects (those whose start is compatible with
the type/maxgap); the exact predicate then filters the candidates, in numpy,
chunked to bound memory.

Predicate (derived from R's behaviour, see NOTES.md). With ``g = maxgap`` and
``m = minoverlap``, ``gg = max(g, 0)``, overlap length
``ol = min(qe, se) - max(qs, ss) + 1``:

* ``any``:    ``qs <= se + g + 1 and ss <= qe + g + 1``  (g = -1: closed-interval overlap)
* ``start``:  ``|qs - ss| <= gg``
* ``end``:    ``|qe - se| <= gg``
* ``equal``:  start and end
* ``within``: ``ss <= qs and qe <= se`` and, if ``g > 0``, ``(se-ss) - (qe-qs) <= g``
              (``g = 0`` means *no* limit, like -1: an R quirk)
* all types additionally require ``ol >= m`` when ``m > 0``.
"""

from __future__ import annotations

import numpy as np

I64 = np.int64
TYPES = ("any", "start", "end", "within", "equal")
SELECTS = ("all", "first", "last", "arbitrary")
_CHUNK = 4_000_000  # max candidate pairs materialised at once


def check_args(type: str, maxgap: int, minoverlap: int) -> None:
    if type not in TYPES:
        raise ValueError(f"'type' must be one of {TYPES}")
    if maxgap < -1:
        raise ValueError("'maxgap' must be >= -1")
    if minoverlap < 0:
        raise ValueError("'minoverlap' cannot be negative")
    if type == "any" and maxgap != -1 and minoverlap != 0:
        raise ValueError(
            "when 'type' is \"any\", at least one of 'maxgap' and 'minoverlap' must be set to its default value"
        )


def predicate(
    qs: np.ndarray, qe: np.ndarray, ss: np.ndarray, se: np.ndarray, type: str, maxgap: int, minoverlap: int
) -> np.ndarray:
    gg = max(maxgap, 0)
    if type == "any":
        ok = (qs <= se + maxgap + 1) & (ss <= qe + maxgap + 1)
    elif type == "start":
        ok = np.abs(qs - ss) <= gg
    elif type == "end":
        ok = np.abs(qe - se) <= gg
    elif type == "equal":
        ok = (np.abs(qs - ss) <= gg) & (np.abs(qe - se) <= gg)
    else:  # within
        ok = (ss <= qs) & (qe <= se)
        if maxgap > 0:
            ok &= (se - ss) - (qe - qs) <= maxgap
    if minoverlap > 0:
        ok &= np.minimum(qe, se) - np.maximum(qs, ss) + 1 >= minoverlap
    return ok


def overlap_pairs(
    qs: np.ndarray,
    qw: np.ndarray,
    ss: np.ndarray,
    sw: np.ndarray,
    type: str = "any",
    maxgap: int = -1,
    minoverlap: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """All ``(query, subject)`` index pairs satisfying the predicate, sorted by
    query then subject."""
    check_args(type, maxgap, minoverlap)
    nq = len(qs)
    if nq == 0 or len(ss) == 0:
        return np.array([], dtype=I64), np.array([], dtype=I64)
    qe = qs + qw - 1
    se = ss + sw - 1
    order = np.argsort(ss, kind="stable")
    ss_sorted = ss[order]
    gg = max(maxgap, 0)
    if type in ("start", "equal"):
        lo = np.searchsorted(ss_sorted, qs - gg, side="left")
        hi = np.searchsorted(ss_sorted, qs + gg, side="right")
    else:
        if type == "any":
            bound = qe + maxgap + 1
        elif type == "end":
            bound = qe + gg + 1  # ss <= se + 1 (empty subjects) <= qe + gg + 1
        else:  # within
            bound = qs
        lo = np.zeros(nq, dtype=I64)
        hi = np.searchsorted(ss_sorted, bound, side="right")
    counts = (hi - lo).astype(I64)
    qi_parts: list[np.ndarray] = []
    si_parts: list[np.ndarray] = []
    # process queries in chunks so that sum(counts) per chunk stays bounded
    cum = np.cumsum(counts)
    begin = 0
    while begin < nq:
        base = cum[begin - 1] if begin else 0
        end = int(np.searchsorted(cum, base + _CHUNK, side="right"))
        end = max(end, begin + 1)
        c = counts[begin:end]
        total = int(c.sum())
        if total:
            q_idx = np.repeat(np.arange(begin, end, dtype=I64), c)
            first = np.cumsum(c) - c
            offs = np.arange(total, dtype=I64) - np.repeat(first, c) + np.repeat(lo[begin:end], c)
            s_idx = order[offs]
            keep = predicate(qs[q_idx], qe[q_idx], ss[s_idx], se[s_idx], type, maxgap, minoverlap)
            qi_parts.append(q_idx[keep])
            si_parts.append(s_idx[keep])
        begin = end
    if not qi_parts:
        return np.array([], dtype=I64), np.array([], dtype=I64)
    qi = np.concatenate(qi_parts)
    si = np.concatenate(si_parts)
    o = np.lexsort((si, qi))
    return qi[o], si[o]


def select_hits(qi: np.ndarray, si: np.ndarray, nq: int, select: str, no_match: int = -1) -> np.ndarray:
    """Reduce sorted ``(qi, si)`` pairs to one subject per query (``select`` first/last/arbitrary).

    ``arbitrary`` returns the first hit (R says any hit may be returned)."""
    if select not in ("first", "last", "arbitrary"):
        raise ValueError(f"'select' must be one of {SELECTS}")
    out = np.full(nq, no_match, dtype=I64)
    if len(qi) == 0:
        return out
    boundary = np.flatnonzero(qi[1:] != qi[:-1])
    if select == "last":
        pos = np.concatenate((boundary, [len(qi) - 1]))
    else:
        pos = np.concatenate(([0], boundary + 1))
    out[qi[pos]] = si[pos]
    return out
