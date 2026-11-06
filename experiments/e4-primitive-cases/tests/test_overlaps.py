"""Overlap kernel vs an O(n*m) pure-Python predicate on seeded random inputs.

This checks the sorted-prefix candidate bound, not the predicate itself;
the predicate is pinned against R by tests/../cases/iranges/overlaps.yaml.
"""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e4 import IRanges


def brute(q: IRanges, s: IRanges, maxgap: int, minoverlap: int, type_: str) -> list[tuple[int, int]]:
    out = []
    for i, (qs, qe) in enumerate(zip(q.start.tolist(), q.end.tolist(), strict=True)):
        for j, (ss, se) in enumerate(zip(s.start.tolist(), s.end.tolist(), strict=True)):
            g = max(maxgap, 0)
            match type_:
                case "any":
                    ok = qs - se - 1 <= maxgap and ss - qe - 1 <= maxgap
                case "start":
                    ok = abs(qs - ss) <= g
                case "end":
                    ok = abs(qe - se) <= g
                case "equal":
                    ok = abs(qs - ss) <= g and abs(qe - se) <= g
                case _:
                    ok = ss <= qs and qe <= se and (maxgap <= 0 or (se - ss) - (qe - qs) <= maxgap)
            if ok and minoverlap > 0:
                ok = min(qe, se) - max(qs, ss) + 1 >= minoverlap
            if ok:
                out.append((i, j))
    return out


def rand(seed: int, n: int) -> IRanges:
    rng = np.random.default_rng(seed)
    return IRanges(start=rng.integers(0, 40, n), width=rng.integers(0, 9, n))


@pytest.mark.parametrize("seed", range(20))
@pytest.mark.parametrize("type_", ["any", "start", "end", "within", "equal"])
@pytest.mark.parametrize("maxgap,minoverlap", [(-1, 0), (0, 0), (3, 0), (-1, 2)])
def test_kernel_matches_brute_force(seed: int, type_: str, maxgap: int, minoverlap: int) -> None:
    if type_ == "any" and maxgap != -1 and minoverlap:
        pytest.skip("invalid combination")
    q, s = rand(seed, 9), rand(seed + 77, 11)
    h = q.find_overlaps(s, maxgap=maxgap, minoverlap=minoverlap, type=type_)
    assert list(zip(h.query.tolist(), h.subject.tolist(), strict=True)) == brute(q, s, maxgap, minoverlap, type_)
    assert (h.query_length, h.subject_length) == (9, 11)


@pytest.mark.parametrize("seed", range(10))
def test_select_consistent_with_all(seed: int) -> None:
    q, s = rand(seed, 9), rand(seed + 5, 11)
    h = q.find_overlaps(s)
    first = q.find_overlaps(s, select="first")
    last = q.find_overlaps(s, select="last")
    arb = q.find_overlaps(s, select="arbitrary")
    for i in range(len(q)):
        hits = h.subject[h.query == i].tolist()
        for sel, pick in ((first, min), (last, max)):
            if hits:
                assert sel[i] == pick(hits)
            else:
                assert sel.mask[i]
        assert (arb[i] in hits) if hits else arb.mask[i]


def test_count_and_any_and_subset_agree() -> None:
    q, s = rand(1, 9), rand(2, 11)
    counts = q.count_overlaps(s)
    assert (q.overlaps_any(s) == (counts > 0)).all()
    assert len(q.subset_by_overlaps(s)) + len(q.subset_by_overlaps(s, invert=True)) == len(q)
