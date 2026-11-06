import warnings

from ranges_e0b import GRanges, IRanges, Seqinfo


def coords(x: GRanges) -> list[tuple[int, int]]:
    return list(zip(x.start.tolist(), x.end.tolist()))


def three() -> GRanges:
    return GRanges("chr1", IRanges(start=[100, 100, 100], end=[110, 110, 110], names=["p", "m", "s"]),
                   strand=["+", "-", "*"], mcols={"v": [1, 2, 3]})


def test_shift_narrow_not_strand_aware() -> None:
    x = three()
    assert coords(x.shift(5)) == [(105, 115)] * 3
    assert coords(x.narrow(start=2, width=3)) == [(101, 103)] * 3
    assert x.shift(5).names == ["p", "m", "s"] and x.shift(5).strand == ["+", "-", "*"]


def test_resize_strand_aware() -> None:
    # R: resize(GRanges("chr1", IRanges(100,110), c("+","-","*")), 3)
    assert coords(three().resize(3)) == [(100, 102), (108, 110), (100, 102)]
    assert coords(three().resize(3, fix="end")) == [(108, 110), (100, 102), (108, 110)]
    assert coords(three().resize(3, fix="center")) == [(104, 106)] * 3


def test_flank_strand_aware() -> None:
    # R: flank(gr, 5)  -- upstream of the 5' end
    assert coords(three().flank(5)) == [(95, 99), (111, 115), (95, 99)]
    assert coords(three().flank(5, start=False)) == [(111, 115), (95, 99), (111, 115)]
    assert coords(three().flank(5, both=True)) == [(95, 104), (106, 115), (95, 104)]
    assert coords(three().flank(-3)) == [(100, 102), (108, 110), (100, 102)]


def test_promoters_strand_aware() -> None:
    # R: promoters(gr, upstream=10, downstream=5)
    assert coords(three().promoters(10, 5)) == [(90, 104), (106, 120), (90, 104)]
    assert three().promoters(10, 5).strand == ["+", "-", "*"]


def test_intra_keeps_seqinfo_mcols() -> None:
    x = three().resize(3)
    assert x.mcols["v"].to_list() == [1, 2, 3] and x.seqlevels == ["chr1"]


def test_flank_out_of_bounds_warns() -> None:
    import pytest

    x = GRanges("chr1", IRanges(start=5, width=10), seqinfo=Seqinfo(["chr1"], [100]))
    with pytest.warns(UserWarning, match="out-of-bound"):
        f = x.flank(10)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        f.trim()


# --- ordering -------------------------------------------------------------------


def mixed() -> GRanges:
    # levels deliberately not alphabetical: chr10 before chr2
    return GRanges(
        ["chr2", "chr10", "chr10", "chr2", "chr10"],
        IRanges(start=[1, 5, 5, 1, 3], width=[2, 2, 2, 2, 2]),
        strand=["+", "-", "+", "*", "+"],
        seqinfo=Seqinfo(["chr10", "chr2"]),
    )


def test_order_seqlevels_then_strand_then_start() -> None:
    x = mixed()
    # chr10: (+3) idx4, (+5) idx2, (-5) idx1 ; chr2: (+1) idx0, (*1) idx3
    assert x.order().tolist() == [4, 2, 1, 0, 3]
    assert x.sort().seqnames == ["chr10", "chr10", "chr10", "chr2", "chr2"]
    assert x.sort().strand == ["+", "+", "-", "+", "*"]


def test_order_follows_seqlevels_not_alphabet() -> None:
    # alphabetical would put chr10 first anyway; reverse the levels to expose the difference
    x = mixed().with_seqinfo(Seqinfo(["chr2", "chr10"]))
    assert x.order().tolist() == [0, 3, 4, 2, 1]


def test_order_decreasing_rank_duplicated_unique() -> None:
    x = GRanges(["a", "a", "b", "a"], IRanges(start=[1, 1, 1, 2], width=1), strand=["+", "+", "+", "+"])
    assert x.duplicated().tolist() == [False, True, False, False]
    assert x.unique().seqnames == ["a", "b", "a"]
    assert x.rank("first").tolist() == [1, 2, 4, 3]
    assert x.order(decreasing=True).tolist() == [2, 3, 0, 1]
    assert x.rank().tolist() == [1.5, 1.5, 4.0, 3.0]


def test_duplicated_strand_matters() -> None:
    x = GRanges("a", IRanges(start=[1, 1], width=1), strand=["+", "-"])
    assert x.duplicated().tolist() == [False, False]
