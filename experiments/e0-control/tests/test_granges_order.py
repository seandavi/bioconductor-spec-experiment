from ranges_e0 import GRanges, IRanges

# seqlevels in order of appearance: chr2, chr1
G = GRanges(
    ["chr2", "chr1", "chr1", "chr1", "chr1"],
    IRanges([1, 5, 5, 5, 3], width=[2, 2, 2, 3, 1]),
    ["+", "-", "-", "-", "*"],
)


def test_order_is_seqlevel_then_strand_then_start_then_width():
    # R: order(g) -> 1 2 3 4 5 : chr2 first (level order, not alphabetical),
    # then on chr1 strands + - * before start
    assert G.order().tolist() == [0, 1, 2, 3, 4]


def test_order_strand_before_start():
    g = GRanges("chr1", IRanges([5, 5, 5, 5], width=2), ["*", "-", "+", "+"])
    # R: order -> 3 4 2 1 ; decreasing -> 1 2 3 4 (ties in original order)
    assert g.order().tolist() == [2, 3, 1, 0]
    assert g.order(decreasing=True).tolist() == [0, 1, 2, 3]
    assert g.sort().strand.tolist() == ["+", "+", "-", "*"]


def test_order_follows_custom_seqlevels():
    from ranges_e0 import Seqinfo

    g = GRanges(["b", "a", "b"], IRanges([5, 3, 1], width=2), seqinfo=Seqinfo(["b", "a", "c"]))
    # R: levels b,a,c -> order 3 1 2
    assert g.order().tolist() == [2, 0, 1]


def test_sort_ignore_strand():
    s = G.sort(ignore_strand=True)
    # R: chr2 1-2 + ; chr1 3 * ; chr1 5-6 - ; 5-6 - ; 5-7 -
    assert s.to_tuples() == [
        ("chr2", 1, 2, "+"), ("chr1", 3, 3, "*"), ("chr1", 5, 6, "-"), ("chr1", 5, 6, "-"), ("chr1", 5, 7, "-")
    ]


def test_sort_decreasing():
    # R: sort(g, decreasing=TRUE): chr1 * 3 ; chr1 - 5-7 ; 5-6 ; 5-6 ; chr2 +
    s = G.sort(decreasing=True)
    assert s.to_tuples() == [
        ("chr1", 3, 3, "*"), ("chr1", 5, 7, "-"), ("chr1", 5, 6, "-"), ("chr1", 5, 6, "-"), ("chr2", 1, 2, "+")
    ]


def test_rank_duplicated_unique():
    # R: rank(ties.method="min") -> 1 2 2 4 5 ; average -> 1 2.5 2.5 4 5
    assert G.rank("min").tolist() == [1, 2, 2, 4, 5]
    assert G.rank("average").tolist() == [1.0, 2.5, 2.5, 4.0, 5.0]
    assert G.rank("first").tolist() == [1, 2, 3, 4, 5]
    assert G.duplicated().tolist() == [False, False, True, False, False]
    assert G.duplicated(from_last=True).tolist() == [False, True, False, False, False]
    assert len(G.unique()) == 4


def test_duplicated_ignores_mcols_and_names_but_not_strand():
    g = GRanges("chr1", IRanges([5, 5], width=2), ["+", "*"], mcols={"s": [1, 1]})
    assert g.duplicated().tolist() == [False, False]  # R: different strand
    assert g.duplicated(ignore_strand=True).tolist() == [False, True]
    h = GRanges("chr1", IRanges([5, 5], width=2), ["+", "+"], mcols={"s": [1, 2]}, names=["a", "b"])
    assert h.duplicated().tolist() == [False, True]


def test_unique_keeps_first_with_mcols_and_names():
    h = GRanges("chr1", IRanges([5, 5, 7], width=2), "+", mcols={"s": [1, 2, 3]}, names=["a", "b", "c"])
    u = h.unique()
    assert u.names == ["a", "c"] and u.mcols["s"].to_list() == [1, 3]


def test_sort_keeps_names_and_seqinfo():
    x = GRanges("chr1", IRanges([5, 1], width=2), "+", names=["a", "b"])
    s = x.sort()
    assert s.names == ["b", "a"] and s.seqlevels == ("chr1",)


def test_empty_ordering():
    e = GRanges.empty()
    assert e.order().tolist() == [] and e.duplicated().tolist() == [] and e.rank().tolist() == []
