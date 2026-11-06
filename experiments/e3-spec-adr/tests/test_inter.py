"""range / reduce / gaps / disjoin (spec/30-inter.md)."""

from __future__ import annotations

import pytest

from helpers import ir, se

Y = ir((1, 3), (4, 6), (8, 8), (9, 12), (20, 25), (5, 5))


@pytest.mark.spec("INTER.1")
def test_range() -> None:
    assert se(ir((5, 4), (1, 0)).range()) == [(1, 4)]
    assert se(ir((1, 3), (5, 4), (10, 13)).range()) == [(1, 13)]
    assert len(ir().range()) == 0


@pytest.mark.spec("INTER.2")
def test_reduce_merges_overlapping_and_adjacent() -> None:
    assert se(Y.reduce()) == [(1, 6), (8, 12), (20, 25)]
    assert se(Y.reduce(min_gapwidth=2)) == [(1, 12), (20, 25)]
    assert se(Y.reduce(min_gapwidth=3)) == [(1, 12), (20, 25)]
    assert se(Y.reduce(min_gapwidth=0)) == [(1, 3), (4, 6), (8, 8), (9, 12), (20, 25)]


@pytest.mark.spec("INTER.2")
def test_reduce_min_gapwidth_zero_keeps_adjacent_apart() -> None:
    assert se(ir((1, 3), (4, 6)).reduce(min_gapwidth=0)) == [(1, 3), (4, 6)]
    assert se(ir((1, 3), (3, 6)).reduce(min_gapwidth=0)) == [(1, 6)]
    assert se(ir((1, 6), (8, 12), (20, 25)).reduce(min_gapwidth=3)) == [(1, 12), (20, 25)]
    with pytest.raises(ValueError):
        ir((1, 3)).reduce(min_gapwidth=-1)


@pytest.mark.spec("INTER.2")
def test_reduce_sorts_unsorted_input_and_empty() -> None:
    assert se(ir((10, 12), (1, 3), (5, 6)).reduce()) == [(1, 3), (5, 6), (10, 12)]
    assert len(ir().reduce()) == 0


@pytest.mark.spec("INTER.3")
def test_reduce_zero_width_ranges() -> None:
    assert se(ir((1, 3), (5, 4), (10, 13), (10, 10)).reduce()) == [(1, 3), (5, 4), (10, 13)]
    assert se(ir((1, 10), (5, 4), (5, 5), (6, 5)).reduce()) == [(1, 10)]
    assert se(ir((1, 3), (4, 3)).reduce()) == [(1, 3)]
    assert se(ir((1, 3), (4, 3)).reduce(min_gapwidth=0)) == [(1, 3), (4, 3)]
    assert se(ir((4, 3), (4, 3), (4, 3)).reduce()) == [(4, 3)]
    assert se(ir((5, 4), (5, 7)).reduce(min_gapwidth=0)) == [(5, 4), (5, 7)]
    assert se(ir((5, 7), (7, 6)).reduce(min_gapwidth=0)) == [(5, 7)]
    assert se(ir((5, 7), (8, 7)).reduce(min_gapwidth=0)) == [(5, 7), (8, 7)]


@pytest.mark.spec("INTER.4")
def test_reduce_drop_empty_ranges() -> None:
    assert se(ir((1, 3), (5, 4), (10, 13)).reduce(drop_empty_ranges=True)) == [(1, 3), (10, 13)]
    assert se(ir((5, 4), (6, 8)).reduce(min_gapwidth=2, drop_empty_ranges=True)) == [(6, 8)]
    assert se(ir((5, 4), (6, 8)).reduce(min_gapwidth=2)) == [(5, 8)]
    assert se(ir((1, 3), (5, 4)).reduce(min_gapwidth=3, drop_empty_ranges=True)) == [(1, 4)]
    assert len(ir((5, 4), (5, 4)).reduce(drop_empty_ranges=True)) == 0


@pytest.mark.spec("INTER.5")
def test_reduce_revmap() -> None:
    r = Y.reduce(with_revmap=True)
    assert r.mcols["revmap"].to_list() == [[0, 1, 5], [2, 3], [4]]
    z = ir((1, 3), (5, 4), (10, 13), (10, 10)).reduce(with_revmap=True)
    assert z.mcols["revmap"].to_list() == [[0], [1], [3, 2]]
    assert ir((10, 12), (1, 3), (5, 6)).reduce(with_revmap=True).mcols["revmap"].to_list() == [[1], [2], [0]]
    ties = ir((3, 6), (3, 5), (3, 6), (3, 5)).reduce(with_revmap=True)
    assert ties.mcols["revmap"].to_list() == [[1, 3, 0, 2]]
    d = ir((10, 12), (5, 4), (1, 3)).reduce(drop_empty_ranges=True, with_revmap=True)
    assert d.mcols["revmap"].to_list() == [[2], [0]]


@pytest.mark.spec("INTER.6")
def test_gaps() -> None:
    x = ir((3, 5), (7, 10), (15, 18))
    assert se(x.gaps()) == [(6, 6), (11, 14)]
    assert se(x.gaps(start=1, end=20)) == [(1, 2), (6, 6), (11, 14), (19, 20)]
    assert se(x.gaps(start=20, end=25)) == [(20, 25)]
    assert se(x.gaps(start=8, end=9)) == []
    assert se(x.gaps(start=12)) == [(12, 14)]
    assert se(x.gaps(end=12)) == [(6, 6), (11, 12)]
    assert se(x.gaps(start=4, end=16)) == [(6, 6), (11, 14)]
    assert se(x.gaps(start=1, end=2)) == [(1, 2)]


@pytest.mark.spec("INTER.6")
def test_gaps_edge_inputs() -> None:
    assert se(ir().gaps(start=1, end=5)) == [(1, 5)]
    assert se(ir().gaps()) == []
    assert se(ir().gaps(start=1)) == []
    assert se(ir((3, 5), (6, 5), (8, 10)).gaps(1, 12)) == [(1, 2), (6, 7), (11, 12)]
    assert se(ir((15, 18), (3, 5)).gaps(1, 20)) == [(1, 2), (6, 14), (19, 20)]
    assert se(ir((1, 5), (3, 8), (10, 12)).gaps(1, 15)) == [(9, 9), (13, 15)]
    assert se(ir((1, 3), (4, 6)).gaps(1, 8)) == [(7, 8)]
    assert se(ir((5, 4)).gaps()) == []
    assert se(ir((5, 4)).gaps(1, 10)) == [(1, 10)]


@pytest.mark.spec("INTER.7")
def test_disjoin() -> None:
    d = ir((1, 5), (3, 8), (10, 12), (5, 4), (20, 19), (6, 6))
    assert se(d.disjoin()) == [(1, 2), (3, 4), (5, 5), (6, 6), (7, 8), (10, 12)]
    assert se(ir((1, 3), (4, 6)).disjoin()) == [(1, 3), (4, 6)]
    assert se(ir((1, 3), (1, 3)).disjoin()) == [(1, 3)]
    assert se(ir((1, 10), (3, 5)).disjoin()) == [(1, 2), (3, 5), (6, 10)]
    assert se(ir((1, 10), (5, 4)).disjoin()) == [(1, 4), (5, 10)]
    assert se(ir((5, 4)).disjoin()) == []
    assert len(ir().disjoin()) == 0
    assert se(ir((6, 8), (1, 3)).disjoin()) == [(1, 3), (6, 8)]


@pytest.mark.spec("INTER.8")
def test_disjoin_revmap() -> None:
    r = ir((1, 5), (3, 8), (10, 12)).disjoin(with_revmap=True)
    assert se(r) == [(1, 2), (3, 5), (6, 8), (10, 12)]
    assert r.mcols["revmap"].to_list() == [[0], [0, 1], [1], [2]]
    z = ir((1, 10), (5, 4), (3, 6)).disjoin(with_revmap=True)
    assert se(z) == [(1, 2), (3, 4), (5, 6), (7, 10)]
    assert z.mcols["revmap"].to_list() == [[0], [0, 2], [0, 2], [0]]
    assert ir((10, 12), (1, 5), (3, 8)).disjoin(with_revmap=True).mcols["revmap"].to_list()[:3] == [[1], [1, 2], [2]]
