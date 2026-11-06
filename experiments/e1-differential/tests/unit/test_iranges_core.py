import numpy as np
import polars as pl
import pytest

from ranges_e1 import IRanges


def test_construct_pairs() -> None:
    a = IRanges(start=[1, 5], end=[3, 4])
    b = IRanges(start=[1, 5], width=[3, 0])
    c = IRanges(end=[3, 4], width=[3, 0])
    for x in (a, b, c):
        assert x.start.tolist() == [1, 5] and x.width.tolist() == [3, 0] and x.end.tolist() == [3, 4]


def test_construct_errors() -> None:
    assert IRanges(start=[1, 2, 3, 7]).width.tolist() == [3, 1]  # R: runs of +1
    with pytest.raises(ValueError):
        IRanges(end=1)
    with pytest.raises(ValueError):
        IRanges(start=5, end=3)
    with pytest.raises(ValueError):
        IRanges(start=1, end=2, width=1)
    with pytest.raises(ValueError):
        IRanges(start=None, end=1, width=None)
    with pytest.raises(ValueError):
        IRanges(start=2**31 - 1, width=2)


def test_immutable_arrays() -> None:
    x = IRanges(start=[1], width=[2])
    with pytest.raises(ValueError):
        x.start[0] = 9


def test_names_padded_and_subset() -> None:
    x = IRanges(start=[1, 2], width=1, names=["a"])
    assert x.names == ("a", None)
    assert x["a"].start.tolist() == [1]
    assert x[[1, 0]].start.tolist() == [2, 1]
    assert x[np.array([True, False])].start.tolist() == [1]
    assert x[-1].start.tolist() == [2]
    with pytest.raises(IndexError):
        x[2]


def test_concat_names_and_mcols() -> None:
    a = IRanges(start=1, width=1, names=["a"], mcols=pl.DataFrame({"x": [1]}))
    b = IRanges(start=2, width=1, mcols=pl.DataFrame({"y": [2.5]}))
    c = IRanges.concat([a, b])
    assert c.names == ("a", "")
    assert c.mcols is not None
    assert c.mcols.columns == ["x", "y"]
    assert c.mcols["x"].to_list() == [1, None]
    assert c.mcols["y"].to_list() == [None, 2.5]
