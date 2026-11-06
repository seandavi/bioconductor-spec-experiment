import numpy as np
import pytest

from ranges_e0 import IRanges

# starts/ends (5,6) (1,5) (1,3) (3,4) (1,3) (5,6)
X = IRanges([5, 1, 1, 3, 1, 5], [6, 5, 3, 4, 3, 6])


def test_order_by_start_then_width():
    # R: order(x) -> 3 5 2 4 1 6   (1-based; here 0-based)
    assert X.order().tolist() == [2, 4, 1, 3, 0, 5]


def test_order_decreasing_keeps_ties_in_original_order():
    # R: order(x, decreasing=TRUE) -> 1 6 4 2 3 5
    assert X.order(decreasing=True).tolist() == [0, 5, 3, 1, 2, 4]


def test_order_zero_width():
    # R: y = IRanges(c(3,3,1), c(2,3,0)) ; order(y) -> 3 1 2
    assert IRanges([3, 3, 1], [2, 3, 0]).order().tolist() == [2, 0, 1]


def test_sort():
    assert X.sort().to_tuples() == [(1, 3), (1, 3), (1, 5), (3, 4), (5, 6), (5, 6)]
    assert X.sort(decreasing=True).to_tuples() == [(5, 6), (5, 6), (3, 4), (1, 5), (1, 3), (1, 3)]
    assert X.sort() == X[X.order()]


def test_sort_carries_names():
    x = X.with_names(list("abcdef"))
    assert x.sort().names == ["c", "e", "b", "d", "a", "f"]


def test_rank_ties():
    # R: rank(x, ties.method=...) -> min: 5 3 1 4 1 5 ; first: 5 3 1 4 2 6 ; max: 6 3 2 4 2 6
    assert X.rank("min").tolist() == [5, 3, 1, 4, 1, 5]
    assert X.rank("first").tolist() == [5, 3, 1, 4, 2, 6]
    assert X.rank("max").tolist() == [6, 3, 2, 4, 2, 6]
    assert X.rank("average").tolist() == [5.5, 3.0, 1.5, 4.0, 1.5, 5.5]
    assert X.rank("last").tolist() == [6, 3, 2, 4, 1, 5]
    assert X.rank().tolist() == X.rank("first").tolist()
    with pytest.raises(ValueError):
        X.rank("random")


def test_duplicated_unique():
    # R: duplicated(x) -> F F F F T T ; fromLast -> T F T F F F
    assert X.duplicated().tolist() == [False, False, False, False, True, True]
    assert X.duplicated(from_last=True).tolist() == [True, False, True, False, False, False]
    assert X.unique().to_tuples() == [(5, 6), (1, 5), (1, 3), (3, 4)]
    assert X.unique(from_last=True).to_tuples() == [(1, 5), (3, 4), (1, 3), (5, 6)]


def test_unique_keeps_first_names():
    x = X.with_names(list("abcdef"))
    assert x.unique().names == ["a", "b", "c", "d"]


def test_duplicated_ignores_names_and_mcols():
    x = IRanges([1, 1], width=2, names=["a", "b"], mcols={"v": [1, 2]})
    assert x.duplicated().tolist() == [False, True]


def test_empty_ordering():
    e = IRanges.empty()
    assert e.order().tolist() == []
    assert len(e.sort()) == 0 and len(e.unique()) == 0
    assert e.rank().tolist() == [] and e.duplicated().tolist() == []


def test_zero_width_ranges_distinct_by_start():
    x = IRanges([5, 5, 6], width=0)
    assert x.duplicated().tolist() == [False, True, False]
