import numpy as np
import polars as pl
import pytest

from ranges_e0 import NO_MATCH, Hits, IRanges, Rle, concat


def test_construct_start_end():
    x = IRanges([1, 5], [3, 4])
    assert x.start.tolist() == [1, 5]
    assert x.end.tolist() == [3, 4]
    assert x.width.tolist() == [3, 0]  # zero-width: end = start - 1
    assert len(x) == 2


def test_construct_start_width_and_end_width():
    assert IRanges(start=[2, 4], width=3).to_tuples() == [(2, 4), (4, 6)]
    assert IRanges(end=[5, 9], width=[2, 0]).to_tuples() == [(4, 5), (10, 9)]


def test_construct_start_only_gives_width_one():
    # R: IRanges(5) -> 5-5
    assert IRanges(5).to_tuples() == [(5, 5)]


def test_construct_all_three_consistent_or_error():
    assert IRanges(1, 5, 5).to_tuples() == [(1, 5)]
    with pytest.raises(ValueError, match="doesn't match"):
        IRanges(1, 5, 3)  # R: "supplied 'width' (3) doesn't match"


def test_construct_errors():
    with pytest.raises(ValueError):
        IRanges(end=5)  # R: at least two of start/end/width
    with pytest.raises(ValueError):
        IRanges(5, 3)  # end < start - 1
    with pytest.raises(ValueError):
        IRanges(1, width=-1)


def test_construct_recycling_and_empty():
    assert IRanges(start=[1, 2, 3], width=2).to_tuples() == [(1, 2), (2, 3), (3, 4)]
    assert len(IRanges.empty()) == 0
    assert len(IRanges([], [])) == 0
    with pytest.raises(ValueError):
        IRanges([1, 2], [5, 6, 7])


def test_construct_integral_float_ok_fractional_rejected():
    assert IRanges(1.0, 3).to_tuples() == [(1, 3)]
    with pytest.raises(ValueError):
        IRanges(1.5, 3)


def test_zero_width_at_start_minus_one():
    x = IRanges(start=[1, 2, 3], width=0)
    assert x.end.tolist() == [0, 1, 2]


def test_names_and_subset_by_name():
    x = IRanges([1, 5, 9], width=2, names=["a", "b", "c"])
    assert x.names == ["a", "b", "c"]
    assert x["b"].to_tuples() == [(5, 6)]
    assert x[["c", "a"]].to_tuples() == [(9, 10), (1, 2)]
    assert x[["c", "a"]].names == ["c", "a"]
    with pytest.raises(KeyError):
        x["z"]  # R: subscript contains invalid names
    with pytest.raises(KeyError):
        IRanges(1, 2)["a"]  # R: cannot subset by character when names are NULL


def test_subset_kinds():
    x = IRanges([1, 5, 9], width=2)
    assert x[0].to_tuples() == [(1, 2)]
    assert x[-1].to_tuples() == [(9, 10)]  # Python semantics, not R's exclusion
    assert x[[2, 0, 0]].to_tuples() == [(9, 10), (1, 2), (1, 2)]
    assert x[np.array([True, False, True])].to_tuples() == [(1, 2), (9, 10)]
    assert x[1:].to_tuples() == [(5, 6), (9, 10)]
    assert x[::-1].to_tuples() == [(9, 10), (5, 6), (1, 2)]
    assert len(x[[]]) == 0
    with pytest.raises(IndexError):
        x[3]
    with pytest.raises(IndexError):
        x[[0, 7]]
    with pytest.raises(IndexError):
        x[np.array([True, False])]


def test_concat_names_and_mcols():
    # R: c(x, IRanges(1,2)) pads names with "" and mcols with NA
    x = IRanges([1, 5], width=[3, 0], names=["a", "b"], mcols={"s": ["u", "v"]})
    y = IRanges(1, 2)
    z = concat(x, y)
    assert z.to_tuples() == [(1, 3), (5, 4), (1, 2)]
    assert z.names == ["a", "b", ""]
    assert z.mcols["s"].to_list() == ["u", "v", None]
    assert concat(IRanges(1, 2), IRanges(3, 4)).names is None
    assert IRanges.concat(y, y).to_tuples() == [(1, 2), (1, 2)]


def test_concat_mcols_aligned_by_name():
    a = IRanges([1], width=1, mcols={"p": [1], "q": ["x"]})
    b = IRanges([2], width=1, mcols={"q": ["y"], "r": [2.5]})
    c = concat(a, b)
    assert c.mcols.columns == ["p", "q", "r"]
    assert c.mcols["p"].to_list() == [1, None]
    assert c.mcols["q"].to_list() == ["x", "y"]
    assert c.mcols["r"].to_list() == [None, 2.5]


def test_mcols_subset_and_with_mcols():
    x = IRanges([1, 5, 9], width=2, mcols={"score": [1.0, 2.0, 3.0]})
    assert x[[2, 0]].mcols["score"].to_list() == [3.0, 1.0]
    assert x.with_mcols(None).mcols.width == 0
    with pytest.raises(ValueError):
        x.with_mcols({"a": [1, 2]})


def test_with_names_and_equality():
    x = IRanges([1, 5], width=2)
    y = x.with_names(["p", "q"])
    assert y.names == ["p", "q"] and x.names is None
    assert x != y
    assert x == IRanges([1, 5], [2, 6])
    assert y.with_names(None) == x


def test_start_is_read_only():
    x = IRanges(1, 5)
    with pytest.raises(ValueError):
        x.start[0] = 9


def test_rle():
    r = Rle.from_array([0, 0, 1, 1, 1, 0])
    assert r.values.tolist() == [0, 1, 0]
    assert r.lengths.tolist() == [2, 3, 1]
    assert len(r) == 6
    assert r.to_numpy().tolist() == [0, 0, 1, 1, 1, 0]
    assert Rle([0, 0, 1], [2, 1, 1]) == Rle([0, 1], [3, 1])  # adjacent equal runs merge
    assert list(Rle([5], [3]).runs()) == [(5, 3)]


def test_hits_sorted_and_validated():
    h = Hits([1, 0, 1], [0, 2, 0 + 1], 2, 3)
    assert h.to_numpy().tolist() == [[0, 2], [1, 0], [1, 1]]
    assert len(h) == 3
    assert list(h) == [(0, 2), (1, 0), (1, 1)]
    with pytest.raises(ValueError):
        Hits([2], [0], 2, 3)


def test_no_match_is_minus_one():
    assert NO_MATCH == -1
