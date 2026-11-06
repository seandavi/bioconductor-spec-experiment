import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from ranges_e0 import NO_MATCH, IRanges

NA = NO_MATCH


def pairs(h):
    return [tuple(p) for p in h.to_numpy().tolist()]


def test_precede_follow_doc_example():
    # R ?precede: query 1-3,3-7,9-10 ; subject 3-3,2-13,10-12
    q = IRanges([1, 3, 9], [3, 7, 10])
    s = IRanges([3, 2, 10], [3, 13, 12])
    assert q.precede(s).tolist() == [2, 2, NA]  # R: 3 3 NA
    assert q.follow(s).tolist() == [NA, NA, 0]  # R: NA NA 1
    assert q.precede().tolist() == [2, 2, NA]  # R: precede(query)
    assert q.follow().tolist() == [NA, NA, 1]  # R: follow(query) -> NA NA 2
    assert IRanges.empty().precede(s).tolist() == []
    assert q.precede(IRanges.empty()).tolist() == [NA] * 3
    assert q.follow(IRanges.empty()).tolist() == [NA] * 3


S = IRanges([10, 10, 5, 20, 5], [12, 15, 6, 25, 5])
Q = IRanges([1, 4, 6, 7, 9, 10, 26], [4, 4, 6, 9, 9, 9, 30])  # incl. a zero-width query (10,9)


def test_precede_ties_and_zero_width():
    # R: precede(Q, S) -> 3 3 1 1 1 1 NA  (ties on start: lowest index)
    assert Q.precede(S).tolist() == [2, 2, 0, 0, 0, 0, NA]
    # R: select="all" -> every subject with the tied start
    assert pairs(Q.precede(S, select="all")) == [
        (0, 2), (0, 4), (1, 2), (1, 4), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (4, 1), (5, 0), (5, 1)
    ]


def test_follow_ties_and_zero_width():
    # R: follow(Q, S) -> NA NA 5 3 3 3 4  (ties on end: highest index)
    assert Q.follow(S).tolist() == [NA, NA, 4, 2, 2, 2, 3]
    assert pairs(Q.follow(S, select="all")) == [(2, 4), (3, 2), (4, 2), (5, 2), (6, 3)]


def test_precede_follow_zero_width_query_and_subject():
    assert IRanges(5, 4).precede(IRanges([5, 6, 3], [7, 8, 4])).tolist() == [0]
    assert IRanges(5, 4).follow(IRanges([5, 6, 3], [7, 8, 4])).tolist() == [2]
    assert IRanges(5, 7).precede(IRanges([8, 8], [7, 9])).tolist() == [0]  # subject (8,7) is empty
    assert IRanges(8, 9).follow(IRanges([5, 4], [7, 3])).tolist() == [0]


def test_precede_follow_self():
    x = IRanges([1, 5, 10], [3, 6, 12])
    assert x.precede().tolist() == [1, 2, NA]
    assert x.follow().tolist() == [NA, 0, 1]


def test_nearest_doc_example():
    # R ?nearest: c(1L, 1L, 3L)
    q = IRanges([1, 3, 9], [2, 7, 10])
    s = IRanges([3, 5, 12], [3, 6, 12])
    assert q.nearest(s).tolist() == [0, 0, 2]


def test_nearest_arbitrary_and_all():
    q = IRanges([1, 3, 9, 5, 15], [2, 7, 10, 6, 16])
    s = IRanges([3, 5, 12], [3, 6, 12])
    # R select="all": (1,1) (2,1) (2,2) (3,3) (4,2) (5,3)  - query 2 overlaps both subjects 1 and 2
    assert pairs(q.nearest(s, select="all")) == [(0, 0), (1, 0), (1, 1), (2, 2), (3, 1), (4, 2)]
    got = q.nearest(s)
    assert got.tolist()[0] == 0 and got.tolist()[2] == 2 and got.tolist()[4] == 2
    assert got.tolist()[1] in (0, 1)  # arbitrary among overlapping


def test_nearest_tie_between_sides():
    # R: IRanges(5,6) vs 1-3 and 8-9 : equal distance -> both with select="all"
    s = IRanges([1, 8], [3, 9])
    assert pairs(IRanges(5, 6).nearest(s, select="all")) == [(0, 0), (0, 1)]
    assert IRanges(5, 6).nearest(s).tolist()[0] in (0, 1)


def test_nearest_adjacent_counts_as_overlap():
    # adjacent (distance 0) wins over a nearer-looking alternative on the other side
    s = IRanges([1, 7], [4, 9])
    assert IRanges(5, 6).nearest(s, select="all").to_numpy().tolist() == [[0, 0], [0, 1]]


def test_nearest_edge_cases():
    assert IRanges(5, 6).nearest(IRanges.empty()).tolist() == [NA]
    assert len(IRanges(5, 6).nearest(IRanges.empty(), select="all")) == 0
    assert IRanges(5, 4).nearest(IRanges([1, 8], [3, 9])).tolist() == [0]  # R: 1
    assert len(IRanges.empty().nearest(IRanges(1, 2))) == 0


def test_nearest_self_excludes_itself():
    x = IRanges([1, 5, 10], [3, 6, 12])
    # R: nearest(x) -> 2 1 2 ; select="all" -> (1,2) (2,1) (3,2)
    assert x.nearest().tolist() == [1, 0, 1]
    assert pairs(x.nearest(select="all")) == [(0, 1), (1, 0), (2, 1)]


def test_distance():
    d = lambda a, b, c, e: int(IRanges(a, b).distance(IRanges(c, e))[0])  # noqa: E731
    # R 4.6 / IRanges 2.47: gap width; 0 when overlapping or adjacent
    assert d(1, 3, 4, 5) == 0
    assert d(1, 3, 5, 6) == 1
    assert d(1, 3, 2, 5) == 0
    assert d(1, 3, 1, 3) == 0
    assert d(5, 6, 1, 3) == 1
    assert d(5, 4, 1, 3) == 1
    assert d(5, 4, 7, 9) == 2
    assert d(5, 4, 5, 7) == 0
    assert d(5, 4, 3, 4) == 0
    assert d(3, 4, 5, 4) == 0
    assert d(5, 4, 3, 10) == 0
    assert d(1, 2, 5, 4) == 2


def test_distance_recycling():
    x = IRanges([1, 10], [2, 12])
    assert x.distance(IRanges(5, 6)).tolist() == [2, 3]  # R: 2 3
    assert IRanges(5, 6).distance(x).tolist() == [2, 3]
    assert len(IRanges.empty().distance(IRanges(5, 6))) == 0
    with pytest.raises(ValueError):
        x.distance(IRanges([1, 2, 3], [4, 5, 6]))


ranges_st = st.lists(st.tuples(st.integers(0, 30), st.integers(0, 6)), min_size=1, max_size=10).map(
    lambda xs: IRanges([a for a, _ in xs], width=[w for _, w in xs])
)


@settings(max_examples=200, deadline=None)
@given(q=ranges_st, s=ranges_st)
def test_precede_follow_nearest_match_naive_definitions(q, s):
    qt, st_ = q.to_tuples(), s.to_tuples()
    for i, (qs, qe) in enumerate(qt):
        after = [(ss, j) for j, (ss, se) in enumerate(st_) if ss > qe]
        exp = min(after)[1] if after else NA  # smallest start, lowest index
        assert q.precede(s).tolist()[i] == exp
        before = [(se, j) for j, (ss, se) in enumerate(st_) if se < qs]
        expf = max(before)[1] if before else NA  # largest end, highest index
        assert q.follow(s).tolist()[i] == expf
        dists = [int(q[i].distance(s[j])[0]) for j in range(len(s))]
        nn = q.nearest(s).tolist()[i]
        assert dists[nn] == min(dists) or min(dists) > 0  # overlapping/adjacent winners have distance 0
    hits = pairs(q.nearest(s, select="all"))
    assert hits == sorted(set(hits))
    assert {i for i, _ in hits} == set(range(len(q)))
