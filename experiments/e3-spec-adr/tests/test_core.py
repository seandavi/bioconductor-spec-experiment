"""IRanges construction, accessors, subsetting, concat (spec/10-iranges-core.md)."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from helpers import ir, se
from ranges_e3 import IRanges, concat


@pytest.mark.spec("CORE.1")
def test_construct_from_any_two() -> None:
    assert se(IRanges(start=[1, 5], end=[3, 4])) == [(1, 3), (5, 4)]
    assert se(IRanges(start=[1, 5], width=[3, 0])) == [(1, 3), (5, 4)]
    assert se(IRanges(end=[3, 4], width=[3, 0])) == [(1, 3), (5, 4)]
    assert IRanges(end=[3, 4], width=[3, 0]).width.tolist() == [3, 0]


@pytest.mark.spec("CORE.1")
@pytest.mark.parametrize("kwargs", [{}, {"start": 1}, {"start": 1, "end": 3, "width": 3}])
def test_exactly_two_required(kwargs: dict[str, int]) -> None:
    with pytest.raises(ValueError):
        IRanges(**kwargs)


@pytest.mark.spec("CORE.2")
@pytest.mark.parametrize(
    "kwargs", [{"start": 3, "end": 1}, {"start": 1, "width": -1}, {"end": 5, "width": -2}]
)
def test_negative_width_refused(kwargs: dict[str, int]) -> None:
    with pytest.raises(ValueError):
        IRanges(**kwargs)


@pytest.mark.spec("CORE.2", "CONV.2")
def test_zero_width_allowed() -> None:
    x = IRanges(start=1, end=0)
    assert x.width.tolist() == [0]


@pytest.mark.spec("CORE.3", "CONV.6")
def test_recycling_and_mismatch() -> None:
    assert se(IRanges(start=[1, 2, 3], width=2)) == [(1, 2), (2, 3), (3, 4)]
    with pytest.raises(ValueError):
        IRanges(start=[1, 2, 3], width=[1, 2])


@pytest.mark.spec("CORE.4")
@pytest.mark.parametrize("bad", [[1, None], [float("nan"), 2.0]])
def test_na_refused(bad: list[object]) -> None:
    with pytest.raises(ValueError):
        IRanges(start=bad, width=1)


@pytest.mark.spec("CORE.5")
def test_empty() -> None:
    x = ir()
    assert len(x) == 0 and x.start.tolist() == [] and se(x[[]]) == []


@pytest.mark.spec("CORE.6", "CONV.1", "CONV.7")
def test_accessors_are_int64() -> None:
    x = IRanges(start=[1, 5, 10], width=[3, 0, 4])
    assert len(x) == 3
    assert x.start.tolist() == [1, 5, 10]
    assert x.end.tolist() == [3, 4, 13]
    assert x.width.tolist() == [3, 0, 4]
    assert x.start.dtype == x.end.dtype == x.width.dtype == np.int64
    big = IRanges(start=2**40, width=2**40)
    assert big.end.tolist() == [2**41 - 1]


@pytest.mark.spec("CORE.7")
def test_names() -> None:
    x = IRanges(start=[1, 2], width=1, names=["a", None])
    assert x.names == ["a", None]
    assert IRanges(start=1, width=1).names is None
    with pytest.raises(ValueError):
        IRanges(start=[1, 2], width=1, names=["a"])


@pytest.mark.spec("CORE.8", "CONV.5")
def test_mcols_and_replacement_is_new_object() -> None:
    x = IRanges(start=[1, 2], width=1)
    assert x.mcols.width == 0
    y = x.with_mcols(pl.DataFrame({"a": [1, 2]}))
    assert y.mcols["a"].to_list() == [1, 2]
    assert x.mcols.width == 0
    with pytest.raises(ValueError):
        x.with_mcols(pl.DataFrame({"a": [1]}))
    z = x.with_names(["p", "q"])
    assert x.names is None and z.names == ["p", "q"]


@pytest.mark.spec("CORE.9", "CONV.3")
def test_subset_by_position_forms() -> None:
    x = ir((1, 3), (5, 4), (10, 13))
    assert se(x[[2, 0]]) == [(10, 13), (1, 3)]
    assert se(x[-1]) == [(10, 13)]
    assert se(x[1]) == [(5, 4)]
    assert se(x[1:]) == [(5, 4), (10, 13)]
    assert se(x[np.array([True, False, True])]) == [(1, 3), (10, 13)]
    assert se(x[[1, 1]]) == [(5, 4), (5, 4)]


@pytest.mark.spec("CORE.9")
@pytest.mark.parametrize("bad", [3, -4, [0, 5], [True, False], "nope"])
def test_subset_refusals(bad: object) -> None:
    x = ir((1, 3), (5, 4), (10, 13), names=["a", "b", "c"])
    with pytest.raises(IndexError):
        x[bad]


@pytest.mark.spec("CORE.9", "CORE.10")
def test_subset_by_name_first_match_and_carries_metadata() -> None:
    x = ir((1, 3), (5, 4), (10, 13), names=["a", "b", "a"])
    x = x.with_mcols(pl.DataFrame({"k": [1, 2, 3]}))
    y = x[["a", "b"]]
    assert se(y) == [(1, 3), (5, 4)]
    assert y.names == ["a", "b"] and y.mcols["k"].to_list() == [1, 2]
    assert x[[2]].names == ["a"]
    assert x["a"].mcols["k"].to_list() == [1]


@pytest.mark.spec("CORE.10")
def test_empty_subset() -> None:
    x = ir((1, 3)).with_mcols(pl.DataFrame({"k": [1]}))
    assert len(x[np.array([], dtype=int)]) == 0
    assert x[[]].mcols.width == 1


@pytest.mark.spec("CORE.11")
def test_concat_order_and_requires_one() -> None:
    assert se(concat(ir((1, 2)), ir((3, 4)), ir((5, 6)))) == [(1, 2), (3, 4), (5, 6)]
    with pytest.raises(ValueError):
        concat()


@pytest.mark.spec("CORE.12")
def test_concat_names() -> None:
    named = ir((1, 2), (3, 4), (5, 6), names=["a", "b", "c"])
    plain = ir((7, 8), (9, 10), (11, 12))
    assert concat(named, plain).names == ["a", "b", "c", "", "", ""]
    assert concat(plain, plain).names is None


@pytest.mark.spec("CORE.12")
def test_concat_mcols_union_with_nulls() -> None:
    a = ir((1, 2)).with_mcols(pl.DataFrame({"x": [1]}))
    b = ir((3, 4)).with_mcols(pl.DataFrame({"y": ["q"]}))
    m = concat(a, b).mcols
    assert m.columns == ["x", "y"]
    assert m["x"].to_list() == [1, None] and m["y"].to_list() == [None, "q"]
    assert concat(a, ir((5, 6))).mcols["x"].to_list() == [1, None]


@pytest.mark.spec("CORE.9")
def test_missing_subscripts_are_refused() -> None:
    x = ir((10, 11), (20, 21))
    with pytest.raises(IndexError):
        x[np.ma.array([1], mask=[True])]
    with pytest.raises(IndexError):
        x[np.ma.array([True, False], mask=[False, True])]
    assert se(x[np.ma.array([1, 0], mask=[False, False])]) == [(20, 21), (10, 11)]
    na_result = ir((1, 3)).find_overlaps(ir((50, 51)), select="first")
    with pytest.raises(IndexError):
        x[na_result]
    mixed = ir((1, 3), (60, 61)).find_overlaps(ir((2, 2)), select="first")
    with pytest.raises(IndexError):
        x[mixed]


@pytest.mark.spec("CORE.9")
def test_boolean_mask_must_match_length_even_when_empty() -> None:
    x = ir((10, 11), (20, 21))
    with pytest.raises(IndexError):
        x[np.array([], dtype=bool)]
    with pytest.raises(IndexError):
        x[np.zeros((2, 1), dtype=bool)]
    assert len(ir()[np.array([], dtype=bool)]) == 0
    assert len(x[np.array([], dtype=int)]) == 0 and len(x[[]]) == 0
