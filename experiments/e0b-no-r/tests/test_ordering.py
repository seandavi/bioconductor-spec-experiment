import numpy as np
from hypothesis import given
from hypothesis import strategies as st

from ranges_e0b import IRanges


def se(x: IRanges) -> list[tuple[int, int]]:
    return list(zip(x.start.tolist(), x.end.tolist()))


X = IRanges(start=[5, 1, 5, 3, 1], end=[9, 4, 6, 3, 4], names=list("abcde"))


def test_order_by_start_then_width() -> None:
    # R: order(IRanges(c(5,1,5,3,1), c(9,4,6,3,4)))  -> 1-based 2 5 4 3 1
    assert X.order().tolist() == [1, 4, 3, 2, 0]


def test_order_decreasing_is_stable_for_ties() -> None:
    # (1,4) appears at 1 and 4; decreasing keeps input order for ties
    assert X.order(decreasing=True).tolist() == [0, 2, 3, 1, 4]


def test_sort_keeps_names() -> None:
    s = X.sort()
    assert se(s) == [(1, 4), (1, 4), (3, 3), (5, 6), (5, 9)]
    assert s.names == ["b", "e", "d", "c", "a"]


def test_rank_methods() -> None:
    # sorted: (1,4)b (1,4)e (3,3)d (5,6)c (5,9)a
    assert X.rank("first").tolist() == [5, 1, 4, 3, 2]
    assert X.rank("min").tolist() == [5, 1, 4, 3, 1]
    assert X.rank("max").tolist() == [5, 2, 4, 3, 2]
    assert X.rank().tolist() == [5.0, 1.5, 4.0, 3.0, 1.5]  # average, UNVERIFIED default


def test_duplicated_unique() -> None:
    # R: duplicated(IRanges(c(1,1,2), c(3,3,3)))
    x = IRanges(start=[1, 1, 2, 1], end=[3, 3, 3, 4])
    assert x.duplicated().tolist() == [False, True, False, False]
    assert se(x.unique()) == [(1, 3), (2, 3), (1, 4)]


def test_empty_ordering() -> None:
    e = IRanges()
    assert e.order().tolist() == [] and e.rank().tolist() == []
    assert e.duplicated().tolist() == [] and len(e.unique()) == 0


@given(st.lists(st.tuples(st.integers(-5, 5), st.integers(0, 4)), max_size=30))
def test_ordering_invariants(pairs: list[tuple[int, int]]) -> None:
    x = IRanges(start=[p[0] for p in pairs], width=[p[1] for p in pairs])
    o = x.order()
    assert sorted(o.tolist()) == list(range(len(x)))
    keys = [(int(x.start[i]), int(x.width[i])) for i in o]
    assert keys == sorted(keys)
    r = x.rank("first")
    assert sorted(r.tolist()) == list(range(1, len(x) + 1))
    assert np.array_equal(r[o], np.arange(1, len(x) + 1))
    seen: set[tuple[int, int]] = set()
    expect = []
    for p in pairs:
        expect.append(p in seen)
        seen.add(p)
    assert x.duplicated().tolist() == expect
