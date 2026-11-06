"""distance / nearest / precede / follow (spec/70-nearest.md)."""

from __future__ import annotations

import pytest

from helpers import ir, na

S = ir((1, 3), (5, 7), (10, 12), (20, 25), (20, 22))
Q = ir((4, 4), (8, 9), (13, 14), (0, 0), (15, 15), (26, 30), (6, 6), (10, 10), (21, 21), (5, 4))


@pytest.mark.spec("NEAR.1")
def test_distance_values() -> None:
    cases = [((1, 3), (4, 6), 0), ((1, 3), (5, 6), 1), ((1, 3), (2, 5), 0), ((5, 4), (1, 3), 1),
             ((11, 10), (1, 10), 0), ((12, 11), (1, 10), 1), ((5, 4), (5, 4), 0), ((5, 4), (6, 5), 1),
             ((5, 4), (9, 8), 4), ((1, 10), (5, 6), 0), ((10, 12), (1, 3), 6)]
    for a, b, d in cases:
        assert ir(a).distance(ir(b)).tolist() == [d], (a, b)
        assert ir(b).distance(ir(a)).tolist() == [d], (b, a)


@pytest.mark.spec("NEAR.1")
def test_distance_recycles_length_one_and_refuses_mismatch() -> None:
    x = ir((1, 3), (1, 3), (10, 12), (5, 4), (5, 6), (1, 10), (8, 7))
    assert x.distance(ir((5, 6))).tolist() == [1, 1, 3, 0, 0, 0, 1]
    assert ir((5, 6)).distance(x).tolist() == [1, 1, 3, 0, 0, 0, 1]
    with pytest.raises(ValueError):
        x.distance(ir((1, 2), (3, 4)))
    assert ir().distance(ir()).tolist() == []


@pytest.mark.spec("NEAR.2")
def test_nearest_all_ties() -> None:
    assert Q.nearest(S, select="all").pairs() == [
        (0, 0), (0, 1), (1, 1), (1, 2), (2, 2), (3, 0), (4, 2), (5, 3), (6, 1), (7, 2), (8, 3), (8, 4), (9, 1),
    ]


@pytest.mark.spec("NEAR.2")
def test_nearest_arbitrary_is_a_tie_member() -> None:
    chosen = na(Q.nearest(S))
    ties = Q.nearest(S, select="all").pairs()
    assert all((q, s) in ties for q, s in enumerate(chosen))
    assert chosen == [0, 1, 2, 0, 2, 3, 1, 2, 3, 1]


@pytest.mark.spec("NEAR.2")
def test_nearest_overlap_and_adjacent_tie_at_zero() -> None:
    assert ir((5, 6)).nearest(ir((6, 8), (3, 4), (8, 9)), select="all").pairs() == [(0, 0), (0, 1)]
    assert ir((5, 6)).nearest(ir((3, 4), (10, 12), (8, 8)), select="all").pairs() == [(0, 0)]
    assert ir((6, 7)).nearest(ir((1, 3), (10, 12), (10, 13), (1, 2)), select="all").pairs() == [(0, 0), (0, 1), (0, 2)]
    assert ir((6, 7)).nearest(ir((1, 3), (9, 12)), select="all").pairs() == [(0, 1)]
    assert ir((5, 5)).nearest(ir((6, 5), (3, 2)), select="all").pairs() == [(0, 0)]
    assert ir((5, 5)).nearest(ir((1, 2), (5, 4)), select="all").pairs() == [(0, 1)]


@pytest.mark.spec("NEAR.2")
def test_nearest_empty_subject_and_self() -> None:
    assert na(ir((1, 3)).nearest(ir())) == [None]
    assert len(ir().nearest(ir((1, 3)), select="all")) == 0
    assert na(S.nearest()) == [1, 0, 1, 4, 3]
    assert na(ir((1, 3)).nearest()) == [None]
    assert na(ir((1, 3), (1, 3), (10, 12)).nearest()) == [1, 0, 0]
    with pytest.raises(ValueError):
        Q.nearest(S, select="first")


@pytest.mark.spec("NEAR.3")
def test_precede() -> None:
    assert na(Q.precede(S)) == [1, 2, 3, 0, 3, None, 2, 3, None, 1]
    assert Q.precede(S, select="all").pairs() == [
        (0, 1), (1, 2), (2, 3), (2, 4), (3, 0), (4, 3), (4, 4), (6, 2), (7, 3), (7, 4), (9, 1),
    ]
    assert na(ir((1, 3)).precede(ir((3, 6)))) == [None]
    assert na(ir((1, 3)).precede(ir((4, 6)))) == [0]
    assert na(ir((1, 3)).precede(ir((4, 6), (5, 6), (3, 6), (4, 9)))) == [0]
    assert ir((1, 3)).precede(ir((4, 6), (5, 6), (3, 6), (4, 9)), select="all").pairs() == [(0, 0), (0, 3)]
    assert na(ir((1, 3)).precede(ir((4, 3), (5, 4)))) == [0]
    assert na(ir((4, 3)).precede(ir((4, 3), (5, 6)))) == [0]
    assert na(S.precede()) == [1, 2, 3, None, None]
    assert na(ir((1, 3)).precede(ir((4, 9), (4, 5), (4, 7)))) == [0]


@pytest.mark.spec("NEAR.4")
def test_follow_first_is_highest_index_of_ties() -> None:
    assert na(Q.follow(S)) == [0, 1, 2, None, 2, 3, 0, 1, 2, 0]
    assert Q.follow(S, select="all").pairs() == [
        (0, 0), (1, 1), (2, 2), (4, 2), (5, 3), (6, 0), (7, 1), (8, 2), (9, 0),
    ]
    x = ir((10, 12))
    assert na(x.follow(ir((1, 3), (3, 5), (7, 9), (8, 9)))) == [3]
    assert x.follow(ir((1, 3), (3, 5), (7, 9), (8, 9)), select="all").pairs() == [(0, 2), (0, 3)]
    assert na(x.follow(ir((8, 9), (7, 9)))) == [1]
    assert na(x.follow(ir((9, 9), (7, 9), (8, 9)))) == [2]
    assert na(x.follow(ir((8, 7), (11, 10)))) == [0]
    assert na(ir((4, 3)).follow(ir((1, 3), (3, 3)))) == [1]
    assert na(ir((10, 12)).follow(ir((1, 30)))) == [None]
    assert na(S.follow()) == [None, 0, 1, 2, 2]
    with pytest.raises(ValueError):
        x.follow(S, select="last")
