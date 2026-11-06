import pytest

from ranges_e0 import GRanges, IRanges, Seqinfo

SI = Seqinfo(["chr1", "chr2", "chr3"], [100, 50, None], [False, False, False], "g1")


def make(seqinfo=SI) -> GRanges:
    # 1 chr2:5-9+ 2 chr1:10-14+ 3 chr1:15-17+ 4 chr1:12-19- 5 chr1:30-34* 6 chr2:8-10+ 7 chr1:40-41-
    return GRanges(
        ["chr2", "chr1", "chr1", "chr1", "chr1", "chr2", "chr1"],
        IRanges([5, 10, 15, 12, 30, 8, 40], width=[5, 5, 3, 8, 5, 3, 2]),
        ["+", "+", "+", "-", "*", "+", "-"],
        seqinfo=seqinfo,
    )


def revmap(g):
    return g.mcols["revmap"].to_list()


def test_reduce_orders_by_seqlevel_then_strand():
    r = make().reduce()
    # R: chr1+ 10-17 ; chr1- 12-19 ; chr1- 40-41 ; chr1* 30-34 ; chr2+ 5-10
    assert r.to_tuples() == [
        ("chr1", 10, 17, "+"), ("chr1", 12, 19, "-"), ("chr1", 40, 41, "-"), ("chr1", 30, 34, "*"), ("chr2", 5, 10, "+")
    ]
    assert r.seqinfo == SI and r.names is None and r.mcols.width == 0


def test_reduce_ignore_strand():
    # R: strand becomes *, groups merge across strands
    r = make().reduce(ignore_strand=True)
    assert r.to_tuples() == [("chr1", 10, 19, "*"), ("chr1", 30, 34, "*"), ("chr1", 40, 41, "*"), ("chr2", 5, 10, "*")]


def test_reduce_min_gapwidth():
    r = make().reduce(min_gapwidth=10)
    assert [t[1:3] for t in r.to_tuples()] == [(10, 17), (12, 19), (40, 41), (30, 34), (5, 10)]
    assert make().reduce(min_gapwidth=21, ignore_strand=True).to_tuples()[0] == ("chr1", 10, 41, "*")


def test_reduce_revmap_indexes_the_input():
    # R (1-based): 2,3 | 4 | 7 | 5 | 1,6
    assert revmap(make().reduce(with_revmap=True)) == [[1, 2], [3], [6], [4], [0, 5]]
    # R ignore.strand: 2,4,3 | 5 | 7 | 1,6  (sweep order by start)
    assert revmap(make().reduce(with_revmap=True, ignore_strand=True)) == [[1, 3, 2], [4], [6], [0, 5]]


def test_reduce_empty_ranges_follow_irange_rules():
    g = GRanges("chr1", IRanges([5, 1], [4, 3]))
    assert g.reduce().to_tuples() == [("chr1", 1, 3, "*"), ("chr1", 5, 4, "*")]
    assert g.reduce(drop_empty_ranges=True).to_tuples() == [("chr1", 1, 3, "*")]


def test_range():
    # R: range(g) -> chr1+ 10-17 ; chr1- 12-41 ; chr1* 30-34 ; chr2+ 5-10
    assert make().range().to_tuples() == [
        ("chr1", 10, 17, "+"), ("chr1", 12, 41, "-"), ("chr1", 30, 34, "*"), ("chr2", 5, 10, "+")
    ]
    assert make().range(ignore_strand=True).to_tuples() == [("chr1", 10, 41, "*"), ("chr2", 5, 10, "*")]


def test_range_with_other_pools_and_merges_seqinfo():
    a = GRanges("chr1", IRanges(1, 5), "+")
    b = GRanges("chr2", IRanges(7, 9), "+")
    with pytest.warns(UserWarning, match="not in the other"):
        r = a.range(b)
    assert r.to_tuples() == [("chr1", 1, 5, "+"), ("chr2", 7, 9, "+")] and r.seqlevels == ("chr1", "chr2")


def test_gaps_known_lengths_cover_every_level_and_strand():
    # R: rows for chr1 & chr2 on all three strands; chr3 (length NA) has none
    assert make().gaps().to_tuples() == [
        ("chr1", 1, 9, "+"), ("chr1", 18, 100, "+"),
        ("chr1", 1, 11, "-"), ("chr1", 20, 39, "-"), ("chr1", 42, 100, "-"),
        ("chr1", 1, 29, "*"), ("chr1", 35, 100, "*"),
        ("chr2", 1, 4, "+"), ("chr2", 11, 50, "+"),
        ("chr2", 1, 50, "-"),
        ("chr2", 1, 50, "*"),
    ]


def test_gaps_unknown_lengths_have_no_trailing_gap_or_empty_groups():
    g = make(Seqinfo(["chr1", "chr2", "chr3"]))
    # R: chr1+ 1-9 ; chr1- 1-11, 20-39 ; chr1* 1-29 ; chr2+ 1-4
    assert g.gaps().to_tuples() == [
        ("chr1", 1, 9, "+"), ("chr1", 1, 11, "-"), ("chr1", 20, 39, "-"), ("chr1", 1, 29, "*"), ("chr2", 1, 4, "+")
    ]


def test_gaps_ignore_strand():
    # R: chr1 1-9, 20-29, 35-39, 42-100 ; chr2 1-4, 11-50 (all strand *)
    assert make().gaps(ignore_strand=True).to_tuples() == [
        ("chr1", 1, 9, "*"), ("chr1", 20, 29, "*"), ("chr1", 35, 39, "*"), ("chr1", 42, 100, "*"),
        ("chr2", 1, 4, "*"), ("chr2", 11, 50, "*"),
    ]


def test_gaps_explicit_start_end_apply_to_every_group():
    # R: gaps(g, start=5, end=20) also reports chr3 (all strands) and empty groups
    assert make().gaps(start=5, end=20).to_tuples() == [
        ("chr1", 5, 9, "+"), ("chr1", 18, 20, "+"), ("chr1", 5, 11, "-"), ("chr1", 20, 20, "-"), ("chr1", 5, 20, "*"),
        ("chr2", 11, 20, "+"), ("chr2", 5, 20, "-"), ("chr2", 5, 20, "*"),
        ("chr3", 5, 20, "+"), ("chr3", 5, 20, "-"), ("chr3", 5, 20, "*"),
    ]


def test_disjoin():
    # R: chr1+ 10-14,15-17 ; chr1- 12-19,40-41 ; chr1* 30-34 ; chr2+ 5-7,8-9,10-10
    assert make().disjoin().to_tuples() == [
        ("chr1", 10, 14, "+"), ("chr1", 15, 17, "+"), ("chr1", 12, 19, "-"), ("chr1", 40, 41, "-"),
        ("chr1", 30, 34, "*"), ("chr2", 5, 7, "+"), ("chr2", 8, 9, "+"), ("chr2", 10, 10, "+"),
    ]
    # R ignore.strand: chr1 10-11,12-14,15-17,18-19,30-34,40-41 ; chr2 5-7,8-9,10-10
    assert [t[1:3] for t in make().disjoin(ignore_strand=True).to_tuples()] == [
        (10, 11), (12, 14), (15, 17), (18, 19), (30, 34), (40, 41), (5, 7), (8, 9), (10, 10)
    ]


def test_disjoin_revmap():
    d = make().disjoin(with_revmap=True)
    assert revmap(d) == [[1], [2], [3], [6], [4], [0], [0, 5], [5]]


def test_empty_object():
    e = GRanges.empty(SI)
    assert len(e.reduce()) == 0 and len(e.disjoin()) == 0 and len(e.range()) == 0
    # empty input still yields whole-sequence gaps where lengths are known
    assert e.gaps().to_tuples()[:3] == [("chr1", 1, 100, "+"), ("chr1", 1, 100, "-"), ("chr1", 1, 100, "*")]


# ---------------------------------------------------------------------------
def pair():
    si = Seqinfo(["chr1", "chr2"], [100, 50], [False, False], "g1")
    x = GRanges(["chr1", "chr1", "chr2", "chr1"], IRanges([1, 20, 5, 50], width=[10, 10, 10, 5]),
                ["+", "-", "+", "*"], seqinfo=si)
    y = GRanges(["chr1", "chr1", "chr2", "chr1"], IRanges([5, 25, 1, 52], width=[10, 10, 3, 10]),
                ["+", "+", "+", "-"], seqinfo=si)
    return x, y


def test_union_is_groupwise_and_star_is_its_own_strand():
    x, y = pair()
    assert x.union(y).to_tuples() == [
        ("chr1", 1, 14, "+"), ("chr1", 25, 34, "+"), ("chr1", 20, 29, "-"), ("chr1", 52, 61, "-"),
        ("chr1", 50, 54, "*"), ("chr2", 1, 3, "+"), ("chr2", 5, 14, "+"),
    ]
    assert x.union(y, ignore_strand=True).to_tuples() == [
        ("chr1", 1, 14, "*"), ("chr1", 20, 34, "*"), ("chr1", 50, 61, "*"), ("chr2", 1, 3, "*"), ("chr2", 5, 14, "*"),
    ]


def test_intersect():
    x, y = pair()
    assert x.intersect(y).to_tuples() == [("chr1", 5, 10, "+")]
    assert x.intersect(y, ignore_strand=True).to_tuples() == [
        ("chr1", 5, 10, "*"), ("chr1", 25, 29, "*"), ("chr1", 52, 54, "*")
    ]


def test_setdiff():
    x, y = pair()
    assert x.setdiff(y).to_tuples() == [
        ("chr1", 1, 4, "+"), ("chr1", 20, 29, "-"), ("chr1", 50, 54, "*"), ("chr2", 5, 14, "+")
    ]
    assert x.setdiff(y, ignore_strand=True).to_tuples() == [
        ("chr1", 1, 4, "*"), ("chr1", 20, 24, "*"), ("chr1", 50, 51, "*"), ("chr2", 5, 14, "*")
    ]


def test_setops_merge_seqinfo_even_for_empty_results():
    x, _ = pair()
    z = GRanges("chrZ", IRanges(1, 5))
    with pytest.warns(UserWarning):
        u = x.union(z)
    assert u.seqlevels == ("chr1", "chr2", "chrZ") and u.seqlengths["chr2"] == 50
    assert u.to_tuples()[-1] == ("chrZ", 1, 5, "*")
    with pytest.warns(UserWarning):
        i = x.intersect(z)
    assert len(i) == 0 and i.seqlevels == ("chr1", "chr2", "chrZ")
    with pytest.warns(UserWarning):
        d = x.setdiff(z)
    assert len(d) == 4 and d.to_tuples()[0] == ("chr1", 1, 10, "+")


def test_setops_conflicting_seqlengths_raise():
    x, _ = pair()
    bad = GRanges("chr1", IRanges(1, 5), seqinfo=Seqinfo(["chr1"], [7]))
    with pytest.raises(ValueError, match="incompatible seqlengths"):
        x.union(bad)


# ---- parallel set ops ----------------------------------------------------------
A = GRanges("chr1", IRanges([1, 5, 10], [4, 8, 12]), ["+", "+", "-"])
B = GRanges("chr1", IRanges([3, 6, 11], [6, 9, 14]), ["+", "*", "-"])


def test_punion_pintersect_psetdiff_star_is_compatible():
    # R: punion -> 1-6+, 5-9+, 10-14- ; pintersect -> 3-4+, 6-8+, 11-12- ; psetdiff -> 1-2+, 5-5+, 10-10-
    assert A.punion(B).to_tuples() == [("chr1", 1, 6, "+"), ("chr1", 5, 9, "+"), ("chr1", 10, 14, "-")]
    assert A.pintersect(B).to_tuples() == [("chr1", 3, 4, "+"), ("chr1", 6, 8, "+"), ("chr1", 11, 12, "-")]
    assert A.psetdiff(B).to_tuples() == [("chr1", 1, 2, "+"), ("chr1", 5, 5, "+"), ("chr1", 10, 10, "-")]


def test_pintersect_adds_hit_column():
    assert A.pintersect(B).mcols["hit"].to_list() == [True, True, True]


def test_punion_incompatible_raises():
    b2 = GRanges("chr1", IRanges([3, 6, 11], [6, 9, 14]), "+")
    with pytest.raises(ValueError, match="compatible"):
        A.punion(b2)
    b3 = GRanges(["chr2", "chr1", "chr1"], IRanges([3, 6, 11], [6, 9, 14]), ["+", "+", "-"])
    with pytest.raises(ValueError, match="compatible"):
        A.punion(b3)
    with pytest.raises(ValueError, match="same length"):
        A.punion(b2[:2])


def test_pintersect_nohit_pairs():
    b2 = GRanges("chr1", IRanges([3, 6, 11], [6, 9, 14]), "+")
    r = A.pintersect(b2)
    # R: third pair has incompatible strands -> empty range at start(x), strand of x, hit=FALSE
    assert r.to_tuples() == [("chr1", 3, 4, "+"), ("chr1", 6, 8, "+"), ("chr1", 10, 9, "-")]
    assert r.mcols["hit"].to_list() == [True, True, False]
    d = A.pintersect(b2, drop_nohit_ranges=True)
    assert d.to_tuples() == [("chr1", 3, 4, "+"), ("chr1", 6, 8, "+")] and d.mcols.width == 0
    b3 = GRanges(["chr2", "chr1", "chr1"], IRanges([3, 6, 11], [6, 9, 14]), ["+", "+", "-"])
    r3 = A.pintersect(b3)
    assert r3.mcols["hit"].to_list() == [False, True, True]
    assert r3.to_tuples()[0] == ("chr1", 1, 0, "+")


def test_pintersect_length_one_y_recycles_and_ignore_strand():
    one = GRanges("chr1", IRanges(3, 9), "*")
    assert [t[1:3] for t in A.pintersect(one).to_tuples()] == [(3, 4), (5, 8), (10, 9)]
    r = A.pintersect(GRanges("chr1", IRanges(3, 11), "+"), ignore_strand=True)
    assert [t[3] for t in r.to_tuples()] == ["+", "+", "-"]  # strand of x
    assert [t[1:3] for t in r.to_tuples()] == [(3, 4), (5, 8), (10, 11)]


def test_psetdiff_incompatible_pairs_keep_x_and_star_resolves_from_y():
    x = GRanges(["chr1", "chr1"], IRanges([1, 1], [10, 10]), ["*", "+"])
    y = GRanges(["chr1", "chr1"], IRanges([1, 1], [5, 5]), ["-", "-"])
    r = x.psetdiff(y)
    assert r.to_tuples() == [("chr1", 6, 10, "-"), ("chr1", 1, 10, "+")]
