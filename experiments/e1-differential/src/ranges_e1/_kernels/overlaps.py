"""All-pairs overlap search on flat arrays. Correctness first: chunked broadcasting."""

from __future__ import annotations

from typing import Any

import numpy as np

from .. import _args

TYPES = ("any", "start", "end", "within", "equal")
_CELLS = 4_000_000  # max query-chunk x subject matrix cells


def normalize_args(type_: Any, select: Any, maxgap: Any, minoverlap: Any) -> tuple[str, str, int, int]:
    """Validate and normalise findOverlaps arguments (partial strings allowed, as match.arg)."""
    t = _args.match_arg(type_, TYPES, "type")
    sel = _args.match_arg(select, ("all", "first", "last", "arbitrary"), "select")
    gap, minov = _args.scalar_int(maxgap, "maxgap"), _args.scalar_int(minoverlap, "minoverlap")
    assert gap is not None and minov is not None
    check_args(t, gap, minov)
    return t, sel, gap, minov


def check_args(type_: str, maxgap: int, minoverlap: int) -> None:
    if type_ not in TYPES:
        raise ValueError(f"'arg' should be one of {', '.join(TYPES)}")
    if maxgap < -1:
        raise ValueError("'maxgap' must be >= -1")
    if minoverlap < 0:
        raise ValueError("'minoverlap' cannot be negative")
    if type_ == "any" and maxgap != -1 and minoverlap != 0:
        raise ValueError("when 'type' is \"any\", at least one of 'maxgap' and 'minoverlap' must be set to its default value")


def uses_query_nclist(nq: int, ns: int) -> bool:
    """R preprocesses the shorter side; this is true when that side is the query."""
    return nq < ns


def find(
    qs: np.ndarray, qw: np.ndarray, ss: np.ndarray, sw: np.ndarray, maxgap: int, minoverlap: int, type_: str,
    ncl_order: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (query_idx, subject_idx) of all matching pairs.

    Sorted by query. Within a query R's order depends on which side it preprocesses
    (an NCList is built on the shorter side, ties -> subject):
    - nq >= ns: subjects in NCList order (start asc, width desc, index asc)
    - nq <  ns: subjects in index order
    `ncl_order` (per subject) overrides the choice: True = NCList order. GRanges decides it
    per seqlevel group, since R runs the engine group by group.

    Predicates were fitted against R on an exhaustive pool incl. empty ranges (see NOTES.md):
    - any:    qs-se-1 <= maxgap and ss-qe-1 <= maxgap
    - start/end/equal: differences <= max(maxgap, 0)
    - within: ss <= qs and qe <= se; for maxgap > 0 also (qs-ss)+(se-qe) <= maxgap
    - minoverlap > 0 additionally needs min(qe,se)-max(qs,ss)+1 >= minoverlap
    """
    nq, ns = len(qs), len(ss)
    empty = np.empty(0, dtype=np.int64)
    if nq == 0 or ns == 0:
        return empty, empty
    if ncl_order is None:
        ncl_order = np.full(ns, nq >= ns)
    perm = np.lexsort((np.arange(ns), np.where(ncl_order, -sw, 0), np.where(ncl_order, ss, 0)))
    ss, sw = ss[perm], sw[perm]
    qe, se = qs + qw - 1, ss + sw - 1
    g = max(maxgap, 0)
    step = max(1, _CELLS // ns)
    out_q: list[np.ndarray] = []
    out_s: list[np.ndarray] = []
    for a in range(0, nq, step):
        b = min(a + step, nq)
        qs_, qe_ = qs[a:b, None], qe[a:b, None]
        if type_ == "any":
            ok = (qs_ - se[None, :] - 1 <= maxgap) & (ss[None, :] - qe_ - 1 <= maxgap)
        elif type_ == "start":
            ok = np.abs(qs_ - ss[None, :]) <= g
        elif type_ == "end":
            ok = np.abs(qe_ - se[None, :]) <= g
        elif type_ == "equal":
            ok = (np.abs(qs_ - ss[None, :]) <= g) & (np.abs(qe_ - se[None, :]) <= g)
        else:
            ok = (ss[None, :] <= qs_) & (qe_ <= se[None, :])
            if maxgap > 0:
                ok &= (qs_ - ss[None, :]) + (se[None, :] - qe_) <= maxgap
        if minoverlap > 0:
            ok &= np.minimum(qe_, se[None, :]) - np.maximum(qs_, ss[None, :]) + 1 >= minoverlap
        qi, si = np.nonzero(ok)
        out_q.append(qi.astype(np.int64) + a)
        out_s.append(perm[si])
    return np.concatenate(out_q), np.concatenate(out_s)
