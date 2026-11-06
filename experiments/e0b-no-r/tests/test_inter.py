from hypothesis import given

from helpers import iranges_st, positions, runs, se
from ranges_e0b import IRanges


def test_range() -> None:
    # R: range(IRanges(c(5,1), c(9,3)))
    assert se(IRanges(start=[5, 1], end=[9, 3]).range()) == [(1, 9)]
    assert len(IRanges().range()) == 0


def test_reduce_basic() -> None:
    # R: reduce(IRanges(c(1,5,20), c(10,12,25)))
    assert se(IRanges(start=[1, 5, 20], end=[10, 12, 25]).reduce()) == [(1, 12), (20, 25)]


def test_reduce_adjacent_merges_by_default() -> None:
    x = IRanges(start=[1, 6], end=[5, 10])
    assert se(x.reduce()) == [(1, 10)]
    assert se(x.reduce(min_gapwidth=0)) == [(1, 5), (6, 10)]


def test_reduce_min_gapwidth() -> None:
    # one uncovered position between the ranges: gap width 1
    x = IRanges(start=[1, 7], end=[5, 10])
    assert se(x.reduce()) == [(1, 5), (7, 10)]
    assert se(x.reduce(min_gapwidth=2)) == [(1, 10)]


def test_reduce_nested_and_unsorted() -> None:
    x = IRanges(start=[10, 1, 3], end=[12, 20, 5])
    assert se(x.reduce()) == [(1, 20)]


def test_reduce_empty_ranges() -> None:
    x = IRanges(start=[5, 1], end=[4, 3])  # empty at 5, and [1,3]
    assert se(x.reduce()) == [(1, 3), (5, 4)]
    assert se(x.reduce(drop_empty_ranges=True)) == [(1, 3)]
    # an empty range inside / adjacent to a real range is absorbed
    y = IRanges(start=[1, 5], end=[10, 4])
    assert se(y.reduce()) == [(1, 10)]


def test_reduce_revmap() -> None:
    x = IRanges(start=[10, 1, 4, 20], end=[15, 5, 8, 22])
    r, rev = x.reduce(with_revmap=True)
    assert se(r) == [(1, 8), (10, 15), (20, 22)]
    assert [m.tolist() for m in rev] == [[1, 2], [0], [3]]


def test_reduce_revmap_with_dropped_empty() -> None:
    x = IRanges(start=[1, 9], end=[3, 8])
    r, rev = x.reduce(drop_empty_ranges=True, with_revmap=True)
    assert se(r) == [(1, 3)] and [m.tolist() for m in rev] == [[0]]


def test_reduce_drops_names_and_mcols() -> None:
    x = IRanges(start=[1], width=2, names=["a"], mcols={"v": [1]})
    r = x.reduce()
    assert r.names is None and r.mcols.width == 0


def test_gaps_default_bounds() -> None:
    # R: gaps(IRanges(c(2,10), c(5,12)))
    x = IRanges(start=[2, 10], end=[5, 12])
    assert se(x.gaps()) == [(6, 9)]


def test_gaps_with_bounds() -> None:
    x = IRanges(start=[2, 10], end=[5, 12])
    # R: gaps(x, start=1, end=15)
    assert se(x.gaps(start=1, end=15)) == [(1, 1), (6, 9), (13, 15)]
    # ranges outside the bounds are clipped away
    assert se(x.gaps(start=4, end=11)) == [(6, 9)]
    assert se(x.gaps(start=20, end=30)) == [(20, 30)]
    assert se(x.gaps(start=1)) == [(1, 1), (6, 9)]
    assert se(x.gaps(end=20)) == [(6, 9), (13, 20)]


def test_gaps_empty_input() -> None:
    assert len(IRanges().gaps()) == 0
    assert se(IRanges().gaps(start=1, end=10)) == [(1, 10)]
    assert len(IRanges().gaps(start=1)) == 0


def test_gaps_overlapping_and_adjacent_inputs() -> None:
    x = IRanges(start=[1, 4, 6], end=[5, 8, 9])
    assert len(x.gaps()) == 0
    assert se(x.gaps(start=0, end=10)) == [(0, 0), (10, 10)]


def test_disjoin() -> None:
    # R: disjoin(IRanges(c(1,5), c(10,12)))
    assert se(IRanges(start=[1, 5], end=[10, 12]).disjoin()) == [(1, 4), (5, 10), (11, 12)]
    assert se(IRanges(start=[1, 4], end=[3, 6]).disjoin()) == [(1, 3), (4, 6)]
    assert se(IRanges(start=[1, 1], end=[3, 3]).disjoin()) == [(1, 3)]
    assert se(IRanges(start=[1, 5], end=[10, 6]).disjoin()) == [(1, 4), (5, 6), (7, 10)]
    assert len(IRanges().disjoin()) == 0


def test_disjoin_gap_not_filled() -> None:
    assert se(IRanges(start=[1, 10], end=[3, 12]).disjoin()) == [(1, 3), (10, 12)]


def test_disjoin_empty_range_splits_but_is_not_output() -> None:
    # UNVERIFIED: empty range adds a breakpoint, as in R's start/end based algorithm
    x = IRanges(start=[1, 5], end=[10, 4])
    assert se(x.disjoin()) == [(1, 4), (5, 10)]


@given(iranges_st())
def test_reduce_invariants(x: IRanges) -> None:
    r = x.reduce(drop_empty_ranges=True)
    assert positions(r) == positions(x)
    assert se(r) == runs(positions(x))  # sorted, disjoint, non-adjacent
    assert r.reduce().equals(r)  # idempotent


@given(iranges_st())
def test_gaps_invariant(x: IRanges) -> None:
    g = x.gaps()
    if len(x) == 0:
        assert len(g) == 0
        return
    lo, hi = int(x.start.min()), int(x.end.max())
    expect = set(range(lo, hi + 1)) - positions(x)
    assert positions(g) == expect
    assert se(g) == runs(expect)


@given(iranges_st())
def test_disjoin_invariant(x: IRanges) -> None:
    d = x.disjoin()
    assert positions(d) == positions(x)
    assert sum(d.width.tolist()) == len(positions(x))  # pieces do not overlap
    assert d.disjoin().equals(d)
    pos = positions(x)
    # each piece is covered identically by the same set of input ranges
    for s, e in se(d):
        cover = {i for i in range(len(x)) if x.start[i] <= s and x.end[i] >= e and x.width[i] > 0}
        assert cover, (s, e)
        for p in range(s, e + 1):
            assert p in pos
            assert {i for i in range(len(x)) if x.width[i] > 0 and x.start[i] <= p <= x.end[i]} == cover


@given(iranges_st())
def test_reduce_revmap_invariant(x: IRanges) -> None:
    r, rev = x.reduce(with_revmap=True)
    assert len(r) == len(rev)
    assert sorted(i for m in rev for i in m.tolist()) == list(range(len(x)))
    for rng_i, m in enumerate(rev):
        assert m.tolist() == sorted(m.tolist())
        assert min(x.start[m]) == r.start[rng_i]
        assert max(x.end[m]) == r.end[rng_i]
