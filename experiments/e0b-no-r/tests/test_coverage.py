import numpy as np
import pytest
from helpers import iranges_st
from hypothesis import given

from ranges_e0b import IRanges, Rle


def test_coverage_basic() -> None:
    # R: coverage(IRanges(c(1,3), c(4,6)))  -> 1 1 2 2 1 1
    c = IRanges(start=[1, 3], end=[4, 6]).coverage()
    assert c.to_numpy().tolist() == [1, 1, 2, 2, 1, 1]
    assert c.values.tolist() == [1, 2, 1] and c.lengths.tolist() == [2, 2, 2]


def test_coverage_gap_and_width() -> None:
    x = IRanges(start=[2, 6], end=[3, 6])
    assert x.coverage().to_numpy().tolist() == [0, 1, 1, 0, 0, 1]
    # R: coverage(x, width=8)
    assert x.coverage(width=8).to_numpy().tolist() == [0, 1, 1, 0, 0, 1, 0, 0]
    # truncating width clips ranges
    assert x.coverage(width=3).to_numpy().tolist() == [0, 1, 1]


def test_coverage_shift_weight() -> None:
    x = IRanges(start=[1, 3], end=[2, 4])
    assert x.coverage(shift=1).to_numpy().tolist() == [0, 1, 1, 1, 1]
    assert x.coverage(weight=[2, 3]).to_numpy().tolist() == [2, 2, 3, 3]
    assert x.coverage(weight=0.5).to_numpy().tolist() == [0.5, 0.5, 0.5, 0.5]


def test_coverage_empty_and_zero_width() -> None:
    assert len(IRanges().coverage()) == 0
    assert len(IRanges(start=5, width=0).coverage()) == 0
    assert IRanges(start=[2, 5], end=[3, 4]).coverage().to_numpy().tolist() == [0, 1, 1]
    assert IRanges().coverage(width=3).to_numpy().tolist() == [0, 0, 0]


def test_coverage_negative_start_clipped() -> None:
    assert IRanges(start=-1, end=2).coverage().to_numpy().tolist() == [1, 1]


def test_coverage_bad_args() -> None:
    with pytest.raises(ValueError):
        IRanges(start=[1, 2], width=1).coverage(weight=[1, 2, 3])
    with pytest.raises(ValueError):
        IRanges(start=1, width=1).coverage(width=-1)


def test_rle_basics() -> None:
    r = Rle.from_array([1, 1, 2, 2, 2, 1])
    assert r.runvalue.tolist() == [1, 2, 1] and r.runlength.tolist() == [2, 3, 1]
    assert len(r) == 6
    assert r == Rle([1, 2, 1], [2, 3, 1])
    assert Rle([1, 1, 2], [1, 1, 3]) == Rle([1, 2], [2, 3])  # equality ignores run splitting
    assert r != Rle([1, 2], [2, 3])
    assert len(Rle.from_array([])) == 0
    with pytest.raises(ValueError):
        Rle([1], [0])
    assert "Rle(length=6" in repr(r)


@given(iranges_st())
def test_coverage_invariants(x: IRanges) -> None:
    c = x.coverage()
    dense = c.to_numpy()
    in_range = sum(max(0, s + w - max(s, 1)) for s, w in zip(x.start.tolist(), x.width.tolist()))
    assert int(dense.sum()) == in_range  # sums to total covered width (positions >= 1)
    expect = np.zeros(len(dense), dtype=np.int64)
    for s, w in zip(x.start.tolist(), x.width.tolist()):
        for p in range(max(s, 1), s + w):
            expect[p - 1] += 1
    assert dense.tolist() == expect.tolist()
    # runs are maximal
    assert not np.any(c.values[1:] == c.values[:-1])
