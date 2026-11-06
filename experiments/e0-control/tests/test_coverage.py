import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from ranges_e0 import IRanges, Rle

X = IRanges([2, 4, 4, 10], [5, 6, 4, 12])


def test_coverage_basic():
    # R: coverage(x) -> lengths 1 2 1 1 1 3 3 ; values 0 1 3 2 1 0 1
    r = X.coverage()
    assert r.lengths.tolist() == [1, 2, 1, 1, 1, 3, 3]
    assert r.values.tolist() == [0, 1, 3, 2, 1, 0, 1]
    assert len(r) == 12
    assert isinstance(r, Rle)


def test_coverage_width_pads_and_truncates():
    # R: width=15 appends 3 zeros ; width=8 cuts the last range
    assert X.coverage(width=15).lengths.tolist() == [1, 2, 1, 1, 1, 3, 3, 3]
    assert X.coverage(width=15).values.tolist() == [0, 1, 3, 2, 1, 0, 1, 0]
    r8 = X.coverage(width=8)
    assert r8.to_numpy().tolist() == [0, 1, 1, 3, 2, 1, 0, 0]
    assert r8.lengths.tolist() == [1, 2, 1, 1, 1, 2]


def test_coverage_shift():
    # R: coverage(x, shift=2) -> length 14, lengths 3 2 1 1 1 3 3
    r = X.coverage(shift=2)
    assert len(r) == 14 and r.lengths.tolist() == [3, 2, 1, 1, 1, 3, 3]
    # R: shift=c(0,1,2,3), width=20
    r = X.coverage(shift=[0, 1, 2, 3], width=20)
    assert r.lengths.tolist() == [1, 3, 2, 1, 5, 3, 5]
    assert r.values.tolist() == [0, 1, 2, 1, 0, 1, 0]


def test_coverage_weight():
    # R: weight=c(1,2,3,4) -> values 0 1 6 3 2 0 4
    r = X.coverage(weight=[1, 2, 3, 4])
    assert r.values.tolist() == [0, 1, 6, 3, 2, 0, 4]
    assert r.lengths.tolist() == [1, 2, 1, 1, 1, 3, 3]
    assert X.coverage(weight=[0.5, 0.5, 0.5, 0.5]).values.dtype.kind == "f"


def test_coverage_empty_ranges_and_inputs():
    # R: coverage(IRanges(5,4)) is 4 zeros (length = max end)
    r = IRanges(5, 4).coverage()
    assert len(r) == 4 and r.values.tolist() == [0]
    assert len(IRanges.empty().coverage()) == 0
    assert len(IRanges.empty().coverage(width=5)) == 5
    assert IRanges([1, 3], [0, 2]).coverage(width=4).values.tolist() == [0]


def test_coverage_clips_at_position_one():
    # R: coverage(IRanges(-2,3)) -> 3 ones ; IRanges(c(-2,3), c(1,5)) -> 1 0 1 1 1
    assert IRanges(-2, 3).coverage().to_numpy().tolist() == [1, 1, 1]
    assert IRanges([-2, 3], [1, 5]).coverage().to_numpy().tolist() == [1, 0, 1, 1, 1]
    assert IRanges(1, 3).coverage(width=2).to_numpy().tolist() == [1, 1]


def test_coverage_runs_are_canonical():
    r = IRanges([1, 3], [2, 4]).coverage()  # adjacent ranges: a single run of 1s
    assert r.lengths.tolist() == [4] and r.values.tolist() == [1]


@settings(max_examples=200, deadline=None)
@given(
    st.lists(st.tuples(st.integers(-3, 30), st.integers(0, 8)), min_size=0, max_size=10),
    st.integers(0, 40) | st.none(),
)
def test_coverage_matches_dense_count(spec, width):
    x = IRanges([a for a, _ in spec], width=[w for _, w in spec]) if spec else IRanges.empty()
    r = x.coverage(width=width)
    length = len(r)
    if width is not None:
        assert length == width
    dense = np.zeros(length, dtype=int)
    for s, e in x.to_tuples():
        for p in range(max(s, 1), min(e, length) + 1):
            dense[p - 1] += 1
    assert r.to_numpy().tolist() == dense.tolist()
    if width is None:
        assert dense.sum() <= int(x.width.sum())
