"""order / sort / rank / duplicated / unique (spec/80-order.md)."""

from __future__ import annotations

import polars as pl
import pytest

from helpers import ir, se

X = ir((5, 9), (1, 4), (5, 8), (3, 3), (1, 4), (5, 9), (5, 4))


@pytest.mark.spec("ORD.1", "ORD.2")
def test_order_ascending_and_decreasing() -> None:
    assert X.order().tolist() == [1, 4, 3, 6, 2, 0, 5]
    assert X.order(decreasing=True).tolist() == [0, 5, 2, 6, 3, 1, 4]
    assert ir((3, 2), (3, 3), (1, 0)).order().tolist() == [2, 0, 1]
    assert ir().order().tolist() == []


@pytest.mark.spec("ORD.1", "ORD.5")
def test_equality_ignores_names_and_mcols() -> None:
    x = ir((1, 3), (1, 3), names=["a", "b"]).with_mcols(pl.DataFrame({"k": [1, 2]}))
    assert x.duplicated().tolist() == [False, True]


@pytest.mark.spec("ORD.3")
def test_sort_carries_metadata() -> None:
    x = ir((3, 4), (1, 2), names=["a", "b"]).with_mcols(pl.DataFrame({"k": ["x", "y"]}))
    s = x.sort()
    assert se(s) == [(1, 2), (3, 4)] and s.names == ["b", "a"] and s.mcols["k"].to_list() == ["y", "x"]
    assert se(X.sort()) == [(1, 4), (1, 4), (3, 3), (5, 4), (5, 8), (5, 9), (5, 9)]
    assert se(X.sort(decreasing=True)) == [(5, 9), (5, 9), (5, 8), (5, 4), (3, 3), (1, 4), (1, 4)]


@pytest.mark.spec("ORD.4")
def test_rank_ties_methods() -> None:
    assert X.rank().tolist() == [6.5, 1.5, 5, 3, 1.5, 6.5, 4]
    assert X.rank("first").tolist() == [6, 1, 5, 3, 2, 7, 4]
    assert X.rank("min").tolist() == [6, 1, 5, 3, 1, 6, 4]
    assert X.rank("max").tolist() == [7, 2, 5, 3, 2, 7, 4]
    assert X.rank("min").dtype.kind == "i" and X.rank().dtype.kind == "f"
    assert ir().rank().tolist() == []
    with pytest.raises(ValueError):
        X.rank("random")


@pytest.mark.spec("ORD.5")
def test_duplicated_and_unique() -> None:
    assert X.duplicated().tolist() == [False, False, False, False, True, True, False]
    assert X.duplicated(from_last=True).tolist() == [True, True, False, False, False, False, False]
    assert se(X.unique()) == [(5, 9), (1, 4), (5, 8), (3, 3), (5, 4)]
    named = ir((1, 3), (1, 3), (2, 2), names=["a", "b", "c"]).unique()
    assert named.names == ["a", "c"]
    assert len(ir().unique()) == 0
