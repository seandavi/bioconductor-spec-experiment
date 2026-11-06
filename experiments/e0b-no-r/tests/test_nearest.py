import pytest
from hypothesis import given, strategies as st

from helpers import iranges_st
from ranges_e0b import NO_HIT, IRanges


def test_distance() -> None:
    # R: distance(IRanges(c(1,10,5), c(3,12,8)), IRanges(c(7,5,6), c(9,11,7)))
    x = IRanges(start=[1, 10, 5, 1], end=[3, 12, 8, 3])
    y = IRanges(start=[7, 5, 6, 4], end=[9, 11, 7, 6])
    assert x.distance(y).tolist() == [3, 0, 0, 0]  # last pair adjacent
    with pytest.raises(ValueError):
        x.distance(y[:2])


def test_precede() -> None:
    # R: precede(IRanges(c(1,10,20),c(5,12,25)), IRanges(c(7,15,30),c(8,16,31)))
    x = IRanges(start=[1, 10, 20], end=[5, 12, 25])
    s = IRanges(start=[7, 15, 30], end=[8, 16, 31])
    assert x.precede(s).tolist() == [0, 1, 2]
    assert IRanges(start=[40], end=[41]).precede(s).tolist() == [NO_HIT]


def test_precede_adjacent_and_overlap() -> None:
    s = IRanges(start=[6, 3], end=[8, 9])
    x = IRanges(start=[1], end=[5])
    # subject 1 starts at 3 (overlaps x) and does not qualify; subject 0 is adjacent and does
    assert x.precede(s).tolist() == [0]


def test_precede_ties_lowest_index_and_all() -> None:
    s = IRanges(start=[10, 10, 20], end=[11, 12, 21])
    x = IRanges(start=[1], end=[5])
    assert x.precede(s).tolist() == [0]
    assert x.precede(s, select="all").pairs() == [(0, 0), (0, 1)]


def test_follow() -> None:
    s = IRanges(start=[1, 3, 10], end=[2, 4, 12])
    x = IRanges(start=[6, 5, 1], end=[7, 5, 3])
    # largest end < start: x0 -> subject 1 (end 4); x1 -> subject 1; x2 -> none
    assert x.follow(s).tolist() == [1, 1, NO_HIT]


def test_follow_ties() -> None:
    s = IRanges(start=[1, 2], end=[4, 4])
    x = IRanges(start=[9], end=[9])
    assert x.follow(s).tolist() == [0]  # UNVERIFIED tie rule: lowest index
    assert x.follow(s, select="all").pairs() == [(0, 0), (0, 1)]


def test_nearest_overlap_wins() -> None:
    s = IRanges(start=[1, 50, 12], end=[5, 60, 20])
    x = IRanges(start=[15, 8, 70], end=[16, 9, 71])
    # x0 overlaps s2; x1: gap to s0 is 2, to s2 is 2 -> tie -> lowest index 0; x2: follows s1
    assert x.nearest(s).tolist() == [2, 0, 1]
    assert x.nearest(s, select="all").pairs() == [(0, 2), (1, 0), (1, 2), (2, 1)]


def test_nearest_adjacent_is_distance_zero() -> None:
    s = IRanges(start=[6, 100], end=[9, 110])
    assert IRanges(start=[1], end=[5]).nearest(s).tolist() == [0]


def test_nearest_self_excludes_self() -> None:
    x = IRanges(start=[1, 10, 14], end=[5, 12, 20])
    assert x.nearest().tolist() == [1, 2, 1]  # 10-12 vs 14-20 gap 1; 1-5 to 10-12 gap 4
    assert x.precede().tolist() == [1, 2, NO_HIT]
    assert x.follow().tolist() == [NO_HIT, 0, 1]


def test_nearest_self_duplicates() -> None:
    x = IRanges(start=[1, 1, 50], end=[5, 5, 60])
    assert x.nearest().tolist() == [1, 0, 0]


def test_nearest_empty_subject() -> None:
    x = IRanges(start=[1], end=[5])
    assert x.nearest(IRanges()).tolist() == [NO_HIT]
    assert x.precede(IRanges()).tolist() == [NO_HIT]
    assert IRanges().nearest(x).tolist() == []


def test_bad_select() -> None:
    x = IRanges(start=[1], end=[5])
    with pytest.raises(ValueError):
        x.nearest(x, select="first")
    with pytest.raises(ValueError):
        x.precede(x, select="arbitrary")


def brute_nearest(q: IRanges, s: IRanges, exclude_self: bool) -> list[list[int]]:
    out = []
    for i in range(len(q)):
        best, ids = None, []
        for j in range(len(s)):
            if exclude_self and i == j:
                continue
            d = max(max(int(q.start[i]), int(s.start[j])) - min(int(q.end[i]), int(s.end[j])) - 1, 0)
            if best is None or d < best:
                best, ids = d, [j]
            elif d == best:
                ids.append(j)
        out.append(ids)
    return out


@given(iranges_st(max_size=8, allow_empty=False), iranges_st(max_size=8, allow_empty=False))
def test_nearest_matches_brute_force(q: IRanges, s: IRanges) -> None:
    got = q.nearest(s, select="all")
    expect = brute_nearest(q, s, False)
    assert got.pairs() == [(i, j) for i, ids in enumerate(expect) for j in ids]


@given(iranges_st(max_size=8, allow_empty=False))
def test_nearest_self_matches_brute_force(q: IRanges) -> None:
    got = q.nearest(select="all")
    expect = brute_nearest(q, q, True)
    assert got.pairs() == [(i, j) for i, ids in enumerate(expect) for j in ids]


@given(iranges_st(max_size=8, allow_empty=False), iranges_st(max_size=8, allow_empty=False))
def test_precede_follow_brute_force(q: IRanges, s: IRanges) -> None:
    pre, fol = q.precede(s), q.follow(s)
    for i in range(len(q)):
        cands = [j for j in range(len(s)) if s.start[j] > q.end[i]]
        exp = min(cands, key=lambda j: (int(s.start[j]), j)) if cands else NO_HIT
        assert pre[i] == exp
        cands = [j for j in range(len(s)) if s.end[j] < q.start[i]]
        exp = min(cands, key=lambda j: (-int(s.end[j]), j)) if cands else NO_HIT
        assert fol[i] == exp
