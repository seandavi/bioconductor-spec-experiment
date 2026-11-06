import pytest
from helpers import iranges_st, positions, runs, se
from hypothesis import given

from ranges_e0b import IRanges


def test_union() -> None:
    # R: union(IRanges(c(1,10), c(5,12)), IRanges(4, 11))
    x = IRanges(start=[1, 10], end=[5, 12])
    assert se(x.union(IRanges(start=4, end=11))) == [(1, 12)]
    assert se(x.union(IRanges(start=20, end=22))) == [(1, 5), (10, 12), (20, 22)]
    assert se(x.union(IRanges())) == [(1, 5), (10, 12)]


def test_union_drops_empty_ranges() -> None:
    assert se(IRanges(start=[1, 9], end=[3, 8]).union(IRanges())) == [(1, 3)]


def test_intersect() -> None:
    # R: intersect(IRanges(1,10), IRanges(5,15))
    assert se(IRanges(start=1, end=10).intersect(IRanges(start=5, end=15))) == [(5, 10)]
    x = IRanges(start=[1, 20], end=[10, 30])
    y = IRanges(start=[5, 25], end=[22, 28])
    assert se(x.intersect(y)) == [(5, 10), (20, 22), (25, 28)]


def test_intersect_adjacent_is_empty_vector() -> None:
    assert len(IRanges(start=1, end=5).intersect(IRanges(start=6, end=10))) == 0


def test_setdiff() -> None:
    # R: setdiff(IRanges(1,10), IRanges(4,6))
    assert se(IRanges(start=1, end=10).setdiff(IRanges(start=4, end=6))) == [(1, 3), (7, 10)]
    assert se(IRanges(start=1, end=10).setdiff(IRanges(start=1, end=10))) == []
    assert se(IRanges(start=1, end=10).setdiff(IRanges())) == [(1, 10)]
    assert se(IRanges().setdiff(IRanges(start=1, end=10))) == []
    assert se(IRanges(start=[1, 3], end=[5, 9]).setdiff(IRanges(start=[2, 8], end=[3, 20]))) == [(1, 1), (4, 7)]


def test_setops_drop_names_and_mcols() -> None:
    x = IRanges(start=[1], width=5, names=["a"], mcols={"v": [1]})
    for r in (x.union(x), x.intersect(x), x.setdiff(IRanges())):
        assert r.names is None and r.mcols.width == 0


def test_punion() -> None:
    # R: punion(IRanges(c(1,10), c(5,12)), IRanges(c(3,11), c(8,20)))
    x = IRanges(start=[1, 10], end=[5, 12])
    y = IRanges(start=[3, 11], end=[8, 20])
    assert se(x.punion(y)) == [(1, 8), (10, 20)]


def test_punion_gap() -> None:
    x, y = IRanges(start=1, end=3), IRanges(start=6, end=8)
    with pytest.raises(ValueError):
        x.punion(y)
    assert se(x.punion(y, fill_gap=True)) == [(1, 8)]
    assert se(x.punion(IRanges(start=4, end=8))) == [(1, 8)]  # adjacent is fine


def test_pintersect() -> None:
    # R: pintersect(IRanges(c(1,10), c(5,20)), IRanges(c(3,12), c(8,15)))
    x = IRanges(start=[1, 10], end=[5, 20], names=["a", "b"])
    y = IRanges(start=[3, 12], end=[8, 15])
    r = x.pintersect(y)
    assert se(r) == [(3, 5), (12, 15)] and r.names == ["a", "b"]


def test_pintersect_non_overlapping() -> None:
    x, y = IRanges(start=1, end=3), IRanges(start=7, end=9)
    with pytest.raises(ValueError):
        x.pintersect(y)
    assert se(x.pintersect(y, resolve_empty="max.start")) == [(7, 6)]
    assert se(x.pintersect(y, resolve_empty="start.x")) == [(1, 0)]
    # UNVERIFIED: adjacent ranges give an empty result without error
    assert se(IRanges(start=1, end=5).pintersect(IRanges(start=6, end=9))) == [(6, 5)]


def test_psetdiff() -> None:
    x = IRanges(start=[1, 1, 1, 5, 1], end=[10, 10, 10, 10, 10])
    y = IRanges(start=[1, 8, 20, 1, 1], end=[4, 12, 30, 12, 10])
    assert se(x.psetdiff(y)) == [(5, 10), (1, 7), (1, 10), (5, 4), (1, 0)]
    with pytest.raises(ValueError):
        IRanges(start=1, end=10).psetdiff(IRanges(start=4, end=6))


def test_parallel_length_mismatch() -> None:
    with pytest.raises(ValueError):
        IRanges(start=[1, 2], width=1).punion(IRanges(start=1, width=1))


@given(iranges_st(), iranges_st())
def test_setop_invariants(x: IRanges, y: IRanges) -> None:
    px, py = positions(x), positions(y)
    for got, expect in (
        (x.union(y), px | py),
        (x.intersect(y), px & py),
        (x.setdiff(y), px - py),
    ):
        assert se(got) == runs(expect)
