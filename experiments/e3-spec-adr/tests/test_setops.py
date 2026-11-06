"""union / intersect / setdiff and parallel versions (spec/40-setops.md)."""

from __future__ import annotations

import pytest

from helpers import ir, se

A = ir((1, 5), (10, 15), (20, 25))
B = ir((4, 11), (12, 13), (30, 35))


@pytest.mark.spec("SET.1")
def test_set_operations_on_covered_positions() -> None:
    assert se(A.union(B)) == [(1, 15), (20, 25), (30, 35)]
    assert se(A.intersect(B)) == [(4, 5), (10, 13)]
    assert se(A.setdiff(B)) == [(1, 3), (14, 15), (20, 25)]
    assert se(B.setdiff(A)) == [(6, 9), (30, 35)]


@pytest.mark.spec("SET.1")
def test_set_operations_adjacent_and_zero_width() -> None:
    assert se(ir((1, 3)).union(ir((4, 6)))) == [(1, 6)]
    assert se(ir((1, 3)).intersect(ir((4, 6)))) == []
    assert se(ir((1, 6)).setdiff(ir((4, 6)))) == [(1, 3)]
    assert se(ir((5, 4)).union(ir((5, 4)))) == []
    assert se(ir((1, 3)).union(ir((4, 3)))) == [(1, 3)]
    assert se(ir((1, 10)).intersect(ir((5, 4), (8, 12)))) == [(8, 10)]
    assert se(ir((1, 10)).setdiff(ir((5, 4)))) == [(1, 10)]
    assert se(ir((5, 4)).setdiff(ir((1, 10)))) == []
    assert se(ir().setdiff(ir((1, 3)))) == []
    assert se(ir((1, 3), (5, 4)).setdiff(ir())) == [(1, 3)]
    assert se(ir((1, 5), (3, 8)).setdiff(ir((2, 2), (7, 7)))) == [(1, 1), (3, 6), (8, 8)]
    assert se(ir((1, 5), (3, 8)).intersect(ir((2, 2), (7, 7)))) == [(2, 2), (7, 7)]


@pytest.mark.spec("SET.2")
def test_parallel_ops_need_equal_length_and_keep_x_names() -> None:
    x = ir((1, 3), names=["a"])
    y = ir((2, 4), names=["b"])
    for op in (x.punion, x.pintersect, x.psetdiff):
        assert op(y).names == ["a"]
        with pytest.raises(ValueError):
            op(ir((1, 2), (3, 4)))


@pytest.mark.spec("SET.3")
def test_punion() -> None:
    assert se(ir((1, 3)).punion(ir((4, 6)))) == [(1, 6)]
    assert se(ir((5, 4)).punion(ir((2, 4)))) == [(2, 4)]
    assert se(ir((5, 4)).punion(ir((5, 8)))) == [(5, 8)]
    assert se(ir((5, 4)).punion(ir((5, 4)))) == [(5, 4)]
    assert se(ir((1, 3)).punion(ir((7, 9)), fill_gap=True)) == [(1, 9)]
    assert se(ir((8, 9)).punion(ir((1, 2)), fill_gap=True)) == [(1, 9)]
    for x, y in ((ir((1, 3)), ir((7, 9))), (ir((5, 4)), ir((1, 3))), (ir((1, 3)), ir((5, 6)))):
        with pytest.raises(ValueError):
            x.punion(y)
    assert se(A.punion(ir((2, 9), (12, 13), (30, 35)), fill_gap=True)) == [(1, 9), (10, 15), (20, 35)]


@pytest.mark.spec("SET.4")
def test_pintersect() -> None:
    assert se(ir((1, 3)).pintersect(ir((4, 6)))) == [(4, 3)]
    assert se(ir((1, 10)).pintersect(ir((5, 4)))) == [(5, 4)]
    assert se(ir((1, 10)).pintersect(ir((11, 10)))) == [(11, 10)]
    assert se(ir((1, 10)).pintersect(ir((1, 0)))) == [(1, 0)]
    assert se(ir((1, 4), (2, 5)).pintersect(ir((2, 9), (2, 9)))) == [(2, 4), (2, 5)]
    x, y = ir((1, 3)), ir((7, 9))
    with pytest.raises(ValueError):
        x.pintersect(y)
    assert se(x.pintersect(y, resolve_empty="max.start")) == [(7, 6)]
    assert se(x.pintersect(y, resolve_empty="start.x")) == [(1, 0)]
    assert se(y.pintersect(x, resolve_empty="start.x")) == [(7, 6)]
    with pytest.raises(ValueError):
        x.pintersect(y, resolve_empty="bogus")


@pytest.mark.spec("SET.5")
def test_psetdiff() -> None:
    assert se(ir((1, 10)).psetdiff(ir((1, 4)))) == [(5, 10)]
    assert se(ir((1, 10)).psetdiff(ir((8, 20)))) == [(1, 7)]
    assert se(ir((1, 3)).psetdiff(ir((7, 9)))) == [(1, 3)]
    assert se(ir((4, 5)).psetdiff(ir((1, 10)))) == [(4, 3)]
    assert se(ir((4, 5)).psetdiff(ir((3, 9)))) == [(4, 3)]
    assert se(ir((4, 8)).psetdiff(ir((5, 20)))) == [(4, 4)]
    assert se(ir((1, 10)).psetdiff(ir((1, 10)))) == [(1, 0)]
    assert se(ir((4, 3)).psetdiff(ir((4, 5)))) == [(4, 3)]


@pytest.mark.spec("SET.5")
def test_psetdiff_refuses_splits() -> None:
    for x, y in (
        (ir((1, 10)), ir((4, 5))),
        (ir((1, 10)), ir((5, 4))),
        (ir((1, 10)), ir((10, 9))),
        (ir((4, 8)), ir((5, 4))),
        (ir((1, 10)), ir((2, 1))),
    ):
        with pytest.raises(ValueError):
            x.psetdiff(y)
    assert se(ir((1, 10)).psetdiff(ir((1, 0)))) == [(1, 10)]
    assert se(ir((1, 10)).psetdiff(ir((11, 10)))) == [(1, 10)]
