import pytest

from ranges_e0 import IRanges


def tup(x):
    return x.to_tuples()


X = IRanges([1, 5, 12], [4, 8, 15])
Y = IRanges([3, 10, 14], [6, 11, 20])


def test_union_intersect_setdiff():
    # R: union(x,y) -> 1-8, 10-20 ; intersect -> 3-6, 14-15 ; setdiff(x,y) -> 1-2, 7-8, 12-13
    assert tup(X.union(Y)) == [(1, 8), (10, 20)]
    assert tup(X.intersect(Y)) == [(3, 6), (14, 15)]
    assert tup(X.setdiff(Y)) == [(1, 2), (7, 8), (12, 13)]
    assert tup(Y.setdiff(X)) == [(10, 11), (16, 20)]


def test_setops_with_empty_ranges():
    # R: empty ranges never appear in the result
    assert tup(IRanges(5, 4).union(IRanges(1, 3))) == [(1, 3)]
    assert len(IRanges(5, 4).union(IRanges(5, 4))) == 0
    assert len(IRanges(1, 10).intersect(IRanges(5, 4))) == 0
    assert tup(IRanges(1, 10).intersect(IRanges(5, 12))) == [(5, 10)]
    assert len(IRanges(1, 4).intersect(IRanges(5, 8))) == 0
    assert tup(IRanges(1, 10).setdiff(IRanges(5, 4))) == [(1, 10)]
    assert len(IRanges(5, 4).setdiff(IRanges(1, 10))) == 0


def test_setdiff_cases():
    assert tup(IRanges(1, 10).setdiff(IRanges(3, 5))) == [(1, 2), (6, 10)]
    assert len(IRanges(1, 10).setdiff(IRanges(1, 10))) == 0
    assert tup(IRanges([1, 20], [10, 30]).setdiff(IRanges(5, 25))) == [(1, 4), (26, 30)]
    assert tup(IRanges([1, 20], [10, 30]).setdiff(IRanges.empty())) == [(1, 10), (20, 30)]
    assert len(IRanges.empty().setdiff(IRanges(5, 25))) == 0


def test_setops_merge_adjacent_and_drop_names():
    x = IRanges([1, 5], [4, 8], names=["a", "b"], mcols={"v": [1, 2]})
    u = x.union(IRanges.empty())
    assert tup(u) == [(1, 8)] and u.names is None and u.mcols.width == 0


A = IRanges([1, 5, 10, 1], [4, 8, 12, 2])
B = IRanges([3, 9, 11, 5], [6, 10, 20, 6])


def test_punion():
    # R: punion(a[1:2], b[3:4], fill.gap=TRUE) -> 1-20, 5-8
    assert tup(A[:2].punion(B[2:], fill_gap=True)) == [(1, 20), (5, 8)]
    assert tup(IRanges(1, 4).punion(IRanges(5, 8))) == [(1, 8)]  # adjacent: no gap
    with pytest.raises(ValueError, match="gap"):
        IRanges(1, 4).punion(IRanges(6, 8))
    assert tup(IRanges(1, 4).punion(IRanges(6, 8), fill_gap=True)) == [(1, 8)]
    with pytest.raises(ValueError, match="same length"):
        IRanges(1, 3).punion(IRanges([1, 2], [4, 4]))  # R: no recycling


def test_punion_names_from_x_else_y():
    x = IRanges(1, 3, names=["a"])
    y = IRanges(2, 5, names=["b"])
    assert x.punion(y).names == ["a"]
    assert IRanges(1, 3).punion(y).names == ["b"]


def test_pintersect():
    assert tup(IRanges(1, 4).pintersect(IRanges(4, 8))) == [(4, 4)]
    assert tup(IRanges(1, 4).pintersect(IRanges(2, 3))) == [(2, 3)]
    # R: adjacent pairs give a definite empty range at the touching position
    assert tup(IRanges(1, 4).pintersect(IRanges(5, 8))) == [(5, 4)]
    assert tup(IRanges(5, 8).pintersect(IRanges(1, 4))) == [(5, 4)]
    # R: disjoint non-adjacent pairs are ambiguous
    with pytest.raises(ValueError, match="ambiguous"):
        IRanges(1, 4).pintersect(IRanges(9, 12))
    assert tup(IRanges(1, 4).pintersect(IRanges(9, 12), resolve_empty="start.x")) == [(1, 0)]
    assert tup(IRanges(1, 4).pintersect(IRanges(9, 12), resolve_empty="max.start")) == [(9, 8)]
    assert tup(A[:2].pintersect(B[2:], resolve_empty="start.x")) == [(1, 0), (5, 6)]


def test_psetdiff():
    # R: psetdiff(a,b) -> 1-2, 5-8, 10-10, 1-2
    assert tup(A.psetdiff(B)) == [(1, 2), (5, 8), (10, 10), (1, 2)]
    assert tup(IRanges(1, 10).psetdiff(IRanges(1, 5))) == [(6, 10)]
    assert tup(IRanges(1, 10).psetdiff(IRanges(5, 10))) == [(1, 4)]
    assert tup(IRanges(1, 10).psetdiff(IRanges(20, 30))) == [(1, 10)]
    assert tup(IRanges(1, 10).psetdiff(IRanges(1, 10))) == [(1, 0)]  # empty at start of x
    assert tup(IRanges(3, 5).psetdiff(IRanges(1, 10))) == [(3, 2)]
    assert tup(IRanges(1, 10).psetdiff(IRanges(5, 12))) == [(1, 4)]
    assert tup(IRanges(5, 10).psetdiff(IRanges(1, 7))) == [(8, 10)]
    with pytest.raises(ValueError, match="strictly inside"):
        IRanges(1, 10).psetdiff(IRanges(3, 5))
