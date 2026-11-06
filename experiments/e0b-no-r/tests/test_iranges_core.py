import numpy as np
import polars as pl
import pytest

from ranges_e0b import IRanges


def test_construct_start_end() -> None:
    # R: IRanges(start=c(1,5), end=c(3,9))
    x = IRanges(start=[1, 5], end=[3, 9])
    assert x.start.tolist() == [1, 5]
    assert x.end.tolist() == [3, 9]
    assert x.width.tolist() == [3, 5]
    assert len(x) == 2


def test_construct_start_width_and_end_width() -> None:
    # R: IRanges(start=c(1,5), width=c(3,5)); IRanges(end=c(3,9), width=c(3,5))
    assert IRanges(start=[1, 5], width=[3, 5]).end.tolist() == [3, 9]
    assert IRanges(end=[3, 9], width=[3, 5]).start.tolist() == [1, 5]


def test_scalar_broadcast() -> None:
    # R: IRanges(start=1:3, width=2)
    x = IRanges(start=[1, 2, 3], width=2)
    assert x.end.tolist() == [2, 3, 4]
    assert len(IRanges(start=5, width=3)) == 1


def test_zero_width() -> None:
    # R: IRanges(start=5, width=0) -> end 4
    x = IRanges(start=5, width=0)
    assert x.end.tolist() == [4]
    assert x.width.tolist() == [0]
    assert len(IRanges(start=[3], end=[2])) == 1


def test_all_three_consistent_and_inconsistent() -> None:
    assert len(IRanges(start=1, end=3, width=3)) == 1
    with pytest.raises(ValueError):
        IRanges(start=1, end=3, width=4)


def test_invalid_construction() -> None:
    with pytest.raises(ValueError):
        IRanges(start=5, end=2)  # negative width
    with pytest.raises(ValueError):
        IRanges(start=[1, 2])  # only one of three
    with pytest.raises(ValueError):
        IRanges(start=[1, 2, 3], end=[4, 5])
    with pytest.raises(ValueError):
        IRanges(start=[1.5], width=[2])


def test_empty() -> None:
    x = IRanges()
    assert len(x) == 0
    assert x.start.tolist() == []


def test_names_and_mcols() -> None:
    x = IRanges(start=[1, 5], width=2, names=["a", "b"], mcols={"score": [10, 20]})
    assert x.names == ["a", "b"]
    assert x.mcols["score"].to_list() == [10, 20]
    with pytest.raises(ValueError):
        IRanges(start=[1, 5], width=2, names=["a"])
    with pytest.raises(ValueError):
        IRanges(start=[1, 5], width=2, mcols={"score": [1, 2, 3]})


def test_immutable_arrays() -> None:
    x = IRanges(start=[1], width=2)
    with pytest.raises(ValueError):
        x.start[0] = 5
    src = np.array([1, 2], dtype=np.int64)
    IRanges(start=src, width=1)
    src[0] = 9  # constructor must not freeze caller's array


def test_subsetting() -> None:
    x = IRanges(start=[1, 5, 10], width=[2, 3, 4], names=["a", "b", "c"], mcols={"v": [1, 2, 3]})
    assert x[1].start.tolist() == [5]
    assert x[-1].names == ["c"]
    assert x[0:2].start.tolist() == [1, 5]
    assert x[[2, 0]].start.tolist() == [10, 1]
    assert x[np.array([True, False, True])].start.tolist() == [1, 10]
    assert x["b"].start.tolist() == [5]
    assert x[["c", "a"]].mcols["v"].to_list() == [3, 1]
    with pytest.raises(IndexError):
        x[3]
    with pytest.raises(ValueError):
        x[np.array([True, False])]
    with pytest.raises(KeyError):
        x["zzz"]


def test_name_lookup_first_match() -> None:
    x = IRanges(start=[1, 5], width=1, names=["a", "a"])
    assert x["a"].start.tolist() == [1]


def test_empty_subset() -> None:
    x = IRanges(start=[1, 5], width=1)
    assert len(x[[]]) == 0
    assert len(x[np.array([False, False])]) == 0


def test_concat() -> None:
    a = IRanges(start=[1], width=2, names=["a"], mcols={"v": [1]})
    b = IRanges(start=[5], width=1, mcols={"v": [2]})
    c = a.concat(b)
    assert c.start.tolist() == [1, 5]
    assert c.names == ["a", ""]
    assert c.mcols["v"].to_list() == [1, 2]
    assert IRanges.combine([a, b, a]).start.tolist() == [1, 5, 1]
    with pytest.raises(ValueError):
        a.concat(IRanges(start=[1], width=1))  # mcols columns differ
    assert len(IRanges.combine([])) == 0


def test_equals() -> None:
    a = IRanges(start=[1, 5], width=2)
    assert a.equals(IRanges(start=[1, 5], end=[2, 6]))
    assert not a.equals(IRanges(start=[1, 5], width=2, names=["x", "y"]))
    assert not a.equals(IRanges(start=[1], width=2))


def test_repr_mentions_length() -> None:
    assert "length 2" in repr(IRanges(start=[1, 5], width=2))
    assert "..." in repr(IRanges(start=np.arange(30), width=1))
