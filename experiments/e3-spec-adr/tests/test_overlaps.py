"""find_overlaps, counts, selects, Hits (spec/60-overlaps.md)."""

from __future__ import annotations

import numpy as np
import pytest

from helpers import ir, na
from ranges_e3.hits import Hits

# subjects: 0:[1,10] 1:[5,4] 2:[5,6] 3:[10,12] 4:[3,3]
S = ir((1, 10), (5, 4), (5, 6), (10, 12), (3, 3))
# queries: 0:[5,4] 1:[5,6] 2:[4,5] 3:[11,11] 4:[0,0] 5:[13,14] 6:[3,2] 7:[1,10]
Q = ir((5, 4), (5, 6), (4, 5), (11, 11), (0, 0), (13, 14), (3, 2), (1, 10))


def pairs(q: object, s: object, **kw: object) -> list[tuple[int, int]]:
    h = q.find_overlaps(s, **kw)  # type: ignore[attr-defined]
    assert isinstance(h, Hits)
    return h.pairs()


@pytest.mark.spec("OVL.1", "OVL.2", "OVL.3")
def test_any_default_including_zero_width() -> None:
    assert pairs(Q, S) == [
        (0, 0), (1, 0), (1, 2), (2, 0), (2, 1), (2, 2), (3, 3),
        (6, 0), (7, 0), (7, 1), (7, 2), (7, 3), (7, 4),
    ]


@pytest.mark.spec("OVL.3")
def test_zero_width_edges_and_pairs() -> None:
    assert pairs(ir((5, 4)), ir((5, 6))) == []
    assert pairs(ir((5, 4)), ir((3, 4))) == []
    assert pairs(ir((5, 4)), ir((1, 10))) == [(0, 0)]
    assert pairs(ir((5, 4)), ir((4, 5))) == [(0, 0)]
    assert pairs(ir((5, 4)), ir((5, 4))) == []


@pytest.mark.spec("OVL.3")
def test_any_max_gap_accepts_adjacent_and_near() -> None:
    assert pairs(ir((5, 8)), ir((1, 3), (2, 4), (10, 12), (11, 12)), max_gap=0) == [(0, 1)]
    assert pairs(ir((5, 8)), ir((1, 3), (2, 4), (10, 12), (11, 12)), max_gap=1) == [(0, 0), (0, 1), (0, 2)]
    assert pairs(ir((5, 8)), ir((1, 3), (2, 4), (10, 12), (11, 13)), max_gap=2) == [(0, 0), (0, 1), (0, 2), (0, 3)]
    assert pairs(ir((1, 5)), ir((7, 9)), max_gap=1) == [(0, 0)]
    assert pairs(ir((1, 5)), ir((6, 9)), max_gap=0) == [(0, 0)]
    assert pairs(ir((5, 4)), ir((5, 6)), max_gap=0) == [(0, 0)]
    assert pairs(ir((5, 4)), ir((5, 4)), max_gap=0) == [(0, 0)]


@pytest.mark.spec("OVL.3", "OVL.6")
def test_any_min_overlap() -> None:
    assert pairs(Q, S, min_overlap=1) == [(1, 0), (1, 2), (2, 0), (2, 2), (3, 3), (7, 0), (7, 2), (7, 3), (7, 4)]
    assert pairs(Q, S, min_overlap=2) == [(1, 0), (1, 2), (2, 0), (7, 0), (7, 2)]
    assert pairs(Q, S, min_overlap=3) == [(7, 0)]
    assert pairs(Q, S, min_overlap=0) == pairs(Q, S)
    assert pairs(ir((5, 8)), ir((1, 5), (4, 5), (6, 7), (8, 12), (9, 12)), min_overlap=2) == [(0, 2)]


@pytest.mark.spec("OVL.2", "OVL.3")
def test_argument_rules() -> None:
    with pytest.raises(ValueError):
        Q.find_overlaps(S, max_gap=0, min_overlap=1)
    with pytest.raises(ValueError):
        Q.find_overlaps(S, max_gap=-2)
    with pytest.raises(ValueError):
        Q.find_overlaps(S, min_overlap=-1)
    with pytest.raises(ValueError):
        Q.find_overlaps(S, type="bogus")
    assert pairs(Q, S, max_gap=-1, min_overlap=2) == pairs(Q, S, min_overlap=2)
    assert pairs(Q, S, type="start", max_gap=0, min_overlap=1) == pairs(Q, S, type="start", min_overlap=1)


@pytest.mark.spec("OVL.4")
def test_start_end_equal() -> None:
    assert pairs(Q, S, type="start") == [(0, 1), (0, 2), (1, 1), (1, 2), (6, 4), (7, 0)]
    assert pairs(Q, S, type="end") == [(0, 1), (1, 2), (7, 0)]
    assert pairs(Q, S, type="equal") == [(0, 1), (1, 2), (7, 0)]
    assert pairs(Q, S, type="start", max_gap=1) == [
        (0, 1), (0, 2), (1, 1), (1, 2), (2, 1), (2, 2), (2, 4), (3, 3), (4, 0), (6, 4), (7, 0),
    ]
    assert pairs(Q, S, type="end", max_gap=1) == [
        (0, 1), (0, 4), (1, 2), (2, 1), (2, 2), (3, 0), (3, 3), (6, 4), (7, 0),
    ]
    sj = ir((5, 8), (5, 9), (6, 8), (1, 20), (4, 8), (3, 9))
    q = ir((5, 8))
    assert pairs(q, sj, type="equal", max_gap=1) == [(0, 0), (0, 1), (0, 2), (0, 4)]
    assert pairs(q, sj, type="equal", max_gap=2) == [(0, 0), (0, 1), (0, 2), (0, 4), (0, 5)]
    assert pairs(q, sj, type="start", max_gap=-1) == pairs(q, sj, type="start", max_gap=0) == [(0, 0), (0, 1)]
    assert pairs(ir((11, 11)), ir((1, 10)), type="start", max_gap=1) == []
    assert pairs(ir((0, 0)), ir((1, 10)), type="start", max_gap=1) == [(0, 0)]


@pytest.mark.spec("OVL.5")
def test_within_and_max_gap_quirk() -> None:
    assert pairs(Q, S, type="within") == [
        (0, 0), (0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (3, 3), (6, 0), (6, 4), (7, 0),
    ]
    subjects = ir((4, 9), (3, 10), (5, 8), (3, 9))
    q = ir((5, 8))
    assert pairs(q, subjects, type="within", max_gap=-1) == [(0, 0), (0, 1), (0, 2), (0, 3)]
    assert pairs(q, subjects, type="within", max_gap=0) == [(0, 0), (0, 1), (0, 2), (0, 3)]
    assert pairs(q, subjects, type="within", max_gap=1) == [(0, 2)]
    assert pairs(q, subjects, type="within", max_gap=2) == [(0, 0), (0, 2)]
    assert pairs(q, subjects, type="within", max_gap=3) == [(0, 0), (0, 2), (0, 3)]
    assert pairs(ir((5, 4)), ir((5, 4), (4, 3)), type="within") == [(0, 0)]


@pytest.mark.spec("OVL.6")
def test_min_overlap_filters_every_type() -> None:
    q = ir((5, 8))
    sj = ir((5, 8), (5, 9), (6, 8), (1, 20), (4, 8), (3, 9))
    assert pairs(q, sj, type="start", min_overlap=4) == [(0, 0), (0, 1)]
    assert pairs(q, sj, type="start", min_overlap=5) == []
    assert pairs(q, sj, type="equal", min_overlap=4) == [(0, 0)]
    assert pairs(q, sj, type="equal", min_overlap=5) == []
    assert pairs(q, ir((5, 4), (9, 12)), type="start", min_overlap=1) == []
    assert pairs(q, ir((5, 4), (9, 12)), type="start") == [(0, 0)]
    assert pairs(ir((5, 4)), ir((4, 5), (5, 6), (3, 8)), type="within", min_overlap=1) == []
    assert pairs(ir((5, 6)), ir((4, 5), (5, 6), (3, 8)), type="within", min_overlap=3) == []


@pytest.mark.spec("OVL.7")
def test_hits_sorted_by_query_then_subject() -> None:
    q = ir((1, 30), (1, 30), (20, 30))
    s = ir((25, 26), (1, 2), (5, 6), (2, 3))
    assert pairs(q, s) == [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3), (2, 0)]


@pytest.mark.spec("OVL.8", "CONV.3")
def test_select_first_last_arbitrary() -> None:
    first = Q.find_overlaps(S, select="first")
    last = Q.find_overlaps(S, select="last")
    arb = Q.find_overlaps(S, select="arbitrary")
    assert isinstance(first, np.ma.MaskedArray)
    assert na(first) == [0, 0, 0, 3, None, None, 0, 0]
    assert na(last) == [0, 2, 2, 3, None, None, 0, 4]
    all_pairs = pairs(Q, S)
    for q, s in enumerate(na(arb)):
        assert (s is None) == all(p[0] != q for p in all_pairs)
        assert s is None or (q, s) in all_pairs
    with pytest.raises(ValueError):
        Q.find_overlaps(S, select="some")


@pytest.mark.spec("OVL.9")
def test_hits_container() -> None:
    h = Q.find_overlaps(S)
    assert isinstance(h, Hits)
    assert len(h) == 13 and h.query_length == 8 and h.subject_length == 5
    assert h.query_hits.dtype == np.int64 and len(h.query_hits) == len(h.subject_hits) == 13
    with pytest.raises(ValueError):
        Hits(np.array([0]), np.array([0, 1]), 1, 2)


@pytest.mark.spec("OVL.10")
def test_count_any_subset() -> None:
    assert Q.count_overlaps(S).tolist() == [1, 2, 3, 1, 0, 0, 1, 5]
    assert Q.count_overlaps(S, min_overlap=2).tolist() == [0, 2, 1, 0, 0, 0, 0, 2]
    assert Q.count_overlaps(S, type="start").tolist() == [2, 2, 0, 0, 0, 0, 1, 1]
    assert Q.overlaps_any(S).tolist() == [True, True, True, True, False, False, True, True]
    sub = Q.subset_by_overlaps(S)
    assert sub.start.tolist() == [5, 5, 4, 11, 3, 1]
    inv = Q.subset_by_overlaps(S, invert=True)
    assert inv.start.tolist() == [0, 13]
    dup = ir((1, 3), (1, 3), (50, 51)).subset_by_overlaps(ir((2, 2)))
    assert dup.start.tolist() == [1, 1]


@pytest.mark.spec("OVL.11")
def test_empty_inputs() -> None:
    assert pairs(ir(), S) == [] and pairs(Q, ir()) == []
    assert na(Q.find_overlaps(ir(), select="first")) == [None] * 8
    assert Q.count_overlaps(ir()).tolist() == [0] * 8
    assert len(ir().find_overlaps(ir())) == 0


@pytest.mark.spec("OVL.2")
def test_large_subject_does_not_hide_hits() -> None:
    s = ir((1, 1000), (500, 501), (900, 900))
    q = ir((950, 960), (499, 499), (1001, 1001))
    assert pairs(q, s) == [(0, 0), (1, 0)]
    assert pairs(q, s, max_gap=0) == [(0, 0), (1, 0), (1, 1), (2, 0)]
