from ranges_e0 import IRanges

# 1-2, 5-6, 20-20, 6-10, 30-31, 31-31
X = IRanges([1, 5, 20, 6, 30, 31], [2, 6, 20, 10, 31, 31])
# with empties: 1-4, (5,4), 5-6, (7,6), (10,9), (3,2)
Z = IRanges([1, 5, 5, 7, 10, 3], [4, 4, 6, 6, 9, 2])


def tup(x):
    return x.to_tuples()


def revmap(x):
    return x.mcols["revmap"].to_list()


def test_reduce_default_merges_overlapping_and_adjacent():
    assert tup(X.reduce()) == [(1, 2), (5, 10), (20, 20), (30, 31)]
    # R: IRanges(c(1,5),c(4,9)) reduce -> 1-9 (adjacent); min.gapwidth=0 -> keeps both
    assert tup(IRanges([1, 5], [4, 9]).reduce()) == [(1, 9)]
    assert tup(IRanges([1, 5], [4, 9]).reduce(min_gapwidth=0)) == [(1, 4), (5, 9)]


def test_reduce_min_gapwidth():
    # R: gap of exactly min.gapwidth is NOT merged
    assert tup(X.reduce(min_gapwidth=3)) == [(1, 10), (20, 20), (30, 31)]
    assert tup(X.reduce(min_gapwidth=10)) == [(1, 31)]
    assert tup(IRanges([1, 10], [3, 12]).reduce(min_gapwidth=6)) == [(1, 3), (10, 12)]  # gap is 6
    assert tup(IRanges([1, 10], [3, 12]).reduce(min_gapwidth=7)) == [(1, 12)]


def test_reduce_drops_names_and_mcols():
    x = IRanges([2, 7], [3, 8], names=["a", "b"], mcols={"m": [1, 2]})
    r = x.reduce()
    assert r.names is None and r.mcols.width == 0


def test_reduce_empty_ranges_default_keeps_isolated_ones():
    # R: reduce(Z) -> 1-6, 10-9  (empties inside/adjacent to real ranges are absorbed)
    assert tup(Z.reduce()) == [(1, 6), (10, 9)]
    assert tup(Z.reduce(drop_empty_ranges=True)) == [(1, 6)]


def test_reduce_empty_ranges_only():
    # R: (5,4),(3,2) stay as two empties; identical empties collapse
    assert tup(IRanges([5, 3], [4, 2]).reduce()) == [(3, 2), (5, 4)]
    assert tup(IRanges([5, 5], [4, 4]).reduce()) == [(5, 4)]
    # empty inside range is absorbed
    assert tup(IRanges([5, 4], [4, 6]).reduce()) == [(4, 6)]


def test_reduce_empty_ranges_respect_min_gapwidth():
    # R: (5,4),(5,5) min.gapwidth=0 -> separate
    assert tup(IRanges([5, 5], [4, 5]).reduce(min_gapwidth=0)) == [(5, 4), (5, 5)]
    # R: (5,4),(6,5): gap 1 -> separate for mgw 0/1, merged into 5-5 for mgw 2
    e = IRanges([5, 6], [4, 5])
    assert tup(e.reduce(min_gapwidth=1)) == [(5, 4), (6, 5)]
    assert tup(e.reduce(min_gapwidth=2)) == [(5, 5)]
    assert tup(e.reduce(min_gapwidth=2, drop_empty_ranges=True)) == []


def test_reduce_revmap():
    # R: revmap 1-based 1 | 2,4 | 3 | 5,6  -> 0-based below
    r = X.reduce(with_revmap=True)
    assert revmap(r) == [[0], [1, 3], [2], [4, 5]]
    assert tup(r) == [(1, 2), (5, 10), (20, 20), (30, 31)]


def test_reduce_revmap_is_in_start_end_order_not_index_order():
    # R: IRanges(c(5,1,3), c(6,2,8)) revmap -> 2 3 1
    assert revmap(IRanges([5, 1, 3], [6, 2, 8]).reduce(with_revmap=True)) == [[1, 2, 0]]
    # ties on start broken by end, then index; R: starts 1,1,1 ends 5,3,4 -> 2 3 1
    assert revmap(IRanges([1, 1, 1], [5, 3, 4]).reduce(with_revmap=True)) == [[1, 2, 0]]
    assert revmap(IRanges([1, 1, 1, 1], [3, 5, 3, 4]).reduce(with_revmap=True)) == [[0, 2, 3, 1]]


def test_reduce_revmap_with_empties():
    # R: two identical empties -> one range, revmap 1,2
    assert revmap(IRanges([5, 5], [4, 4]).reduce(with_revmap=True)) == [[0, 1]]
    r = IRanges([5, 1], [4, 3]).reduce(with_revmap=True, drop_empty_ranges=True)
    assert tup(r) == [(1, 3)] and revmap(r) == [[1]]


def test_reduce_empty_input():
    r = IRanges.empty().reduce(with_revmap=True)
    assert len(r) == 0 and revmap(r) == []


def test_range():
    assert tup(X.range()) == [(1, 31)]
    assert tup(IRanges([5, 3], [4, 2]).range()) == [(3, 4)]
    assert len(IRanges.empty().range()) == 0
    # R: range(x, y) pools both arguments
    assert tup(IRanges([1, 5], [2, 6]).range(IRanges(10, 12))) == [(1, 12)]
    assert IRanges([2, 7], [3, 8], names=["a", "b"]).range().names is None


def test_gaps_default_is_span_of_input():
    assert tup(X.gaps()) == [(3, 4), (11, 19), (21, 29)]
    assert tup(IRanges([2, 7], [3, 8]).gaps()) == [(4, 6)]


def test_gaps_start_end():
    assert tup(X.gaps(start=0)) == [(0, 0), (3, 4), (11, 19), (21, 29)]
    assert tup(X.gaps(end=40)) == [(3, 4), (11, 19), (21, 29), (32, 40)]
    assert tup(X.gaps(start=0, end=40)) == [(0, 0), (3, 4), (11, 19), (21, 29), (32, 40)]
    assert tup(X.gaps(start=5, end=25)) == [(11, 19), (21, 25)]
    assert tup(X.gaps(start=10, end=15)) == [(11, 15)]
    assert tup(IRanges([2, 7], [3, 8]).gaps(start=5)) == [(5, 6)]
    assert tup(IRanges([2, 7], [3, 8]).gaps(end=6)) == [(4, 6)]


def test_gaps_edge_cases():
    assert tup(IRanges.empty().gaps(start=1, end=5)) == [(1, 5)]
    assert len(IRanges.empty().gaps()) == 0
    assert len(X.gaps(start=2, end=1)) == 0
    assert tup(IRanges(10, 12).gaps(start=1, end=5)) == [(1, 5)]  # range beyond the window
    assert len(IRanges(5, 4).gaps()) == 0
    assert len(IRanges([1, 3], [5, 4]).gaps(start=3, end=3)) == 0
    assert tup(IRanges([1, 5], [4, 4]).gaps(start=1, end=10)) == [(5, 10)]  # empty range ignored
    assert len(IRanges([5, 10], [4, 9]).gaps()) == 0


def test_disjoin():
    assert tup(X.disjoin()) == [(1, 2), (5, 5), (6, 6), (7, 10), (20, 20), (30, 30), (31, 31)]
    assert tup(IRanges([1, 3, 3], [10, 5, 3]).disjoin()) == [(1, 2), (3, 3), (4, 5), (6, 10)]
    assert tup(IRanges([1, 1], [3, 3]).disjoin()) == [(1, 3)]
    assert tup(IRanges([1, 4], [3, 6]).disjoin()) == [(1, 3), (4, 6)]  # adjacent stay apart


def test_disjoin_empty_ranges_split_but_are_not_covered():
    # R: disjoin(Z) -> 1-2, 3-4, 5-6  ((3,2) splits 1-4; (10,9) vanishes)
    assert tup(Z.disjoin()) == [(1, 2), (3, 4), (5, 6)]
    assert len(IRanges([5, 5], [4, 4]).disjoin()) == 0
    assert len(IRanges.empty().disjoin()) == 0


def test_disjoin_revmap_ascending_index():
    d = X.disjoin(with_revmap=True)
    assert revmap(d) == [[0], [1], [1, 3], [3], [2], [4], [4, 5]]
    # R: IRanges(c(6,1,3), c(9,10,5)) -> 2 | 2 3 | 1 2 | 2  (1-based)
    d2 = IRanges([6, 1, 3], [9, 10, 5]).disjoin(with_revmap=True)
    assert tup(d2) == [(1, 2), (3, 5), (6, 9), (10, 10)]
    assert revmap(d2) == [[1], [1, 2], [0, 1], [1]]
    d3 = IRanges([1, 4], [3, 6]).disjoin(with_revmap=True)
    assert revmap(d3) == [[0], [1]]


def test_gaps_default_bounds_come_from_non_empty_ranges():
    # R: gaps(IRanges(c(0,3,9,20), c(-1,5,8,25))) -> 6-19 (empties at 0 and 9 don't set bounds)
    x = IRanges([0, 3, 9, 20], [-1, 5, 8, 25])
    assert tup(x.gaps()) == [(6, 19)]
    assert tup(x.gaps(start=1)) == [(1, 2), (6, 19)]
    # R: only an empty range: a missing bound gives nothing, both given -> the whole window
    assert len(IRanges(10, 9).gaps(start=1)) == 0
    assert len(IRanges(10, 9).gaps(end=20)) == 0
    assert tup(IRanges(10, 9).gaps(start=1, end=20)) == [(1, 20)]
