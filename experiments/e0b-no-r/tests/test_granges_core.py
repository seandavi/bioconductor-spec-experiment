import warnings

import pytest

from ranges_e0b import GRanges, IRanges, Seqinfo


def gr() -> GRanges:
    return GRanges(
        ["chr2", "chr1", "chr2"],
        IRanges(start=[1, 5, 10], width=[3, 4, 5], names=["a", "b", "c"]),
        strand=["+", "-", "*"],
        mcols={"score": [1, 2, 3]},
    )


def test_construct_accessors() -> None:
    # R: GRanges(c("chr2","chr1","chr2"), IRanges(c(1,5,10), width=c(3,4,5)), c("+","-","*"))
    x = gr()
    assert x.seqnames == ["chr2", "chr1", "chr2"]
    assert x.strand == ["+", "-", "*"]
    assert x.start.tolist() == [1, 5, 10] and x.end.tolist() == [3, 8, 14] and x.width.tolist() == [3, 4, 5]
    assert x.names == ["a", "b", "c"]
    assert x.mcols["score"].to_list() == [1, 2, 3]
    assert len(x) == 3 and len(x.ranges) == 3
    assert x.ranges.mcols.width == 0


def test_default_seqlevels_are_first_appearance() -> None:
    assert gr().seqlevels == ["chr2", "chr1"]


def test_scalar_seqname_and_strand() -> None:
    x = GRanges("chr1", IRanges(start=[1, 2], width=1))
    assert x.seqnames == ["chr1", "chr1"] and x.strand == ["*", "*"]
    assert GRanges("chr1", IRanges(start=[1, 2], width=1), "-").strand == ["-", "-"]


def test_invalid_inputs() -> None:
    r = IRanges(start=[1, 2], width=1)
    with pytest.raises(ValueError):
        GRanges(["a"], r)
    with pytest.raises(ValueError):
        GRanges("a", r, strand="?")
    with pytest.raises(ValueError):
        GRanges("a", r, strand=["+"] * 3)
    with pytest.raises(ValueError):
        GRanges("a", r, seqinfo=Seqinfo(["b"]))
    with pytest.raises(ValueError):
        GRanges("a", r, mcols={"v": [1]})


def test_seqinfo_accessors() -> None:
    si = Seqinfo(["chr1", "chr2", "chrM"], [100, 200, 16], [False, False, True], "hg38")
    x = GRanges(["chr1"], IRanges(start=1, width=5), seqinfo=si)
    assert x.seqlevels == ["chr1", "chr2", "chrM"]
    assert x.seqlengths == {"chr1": 100, "chr2": 200, "chrM": 16}
    assert x.is_circular == {"chr1": False, "chr2": False, "chrM": True}
    assert x.genome == {"chr1": "hg38", "chr2": "hg38", "chrM": "hg38"}
    assert x.seqinfo is si


def test_subsetting() -> None:
    x = gr()
    y = x[[2, 0]]
    assert y.seqnames == ["chr2", "chr2"] and y.strand == ["*", "+"] and y.names == ["c", "a"]
    assert y.mcols["score"].to_list() == [3, 1]
    assert y.seqlevels == ["chr2", "chr1"]  # seqinfo is preserved even when a level is unused
    assert x["b"].seqnames == ["chr1"]
    assert x[-1].start.tolist() == [10]
    assert x[x.start > 2].names == ["b", "c"]


def test_concat_merges_seqinfo() -> None:
    a = GRanges("chr1", IRanges(start=1, width=5), seqinfo=Seqinfo(["chr1", "chr2"], [100, None]))
    b = GRanges("chr2", IRanges(start=1, width=5), seqinfo=Seqinfo(["chr2", "chr3"], [200, 300]))
    c = a.concat(b)
    assert c.seqnames == ["chr1", "chr2"]
    assert c.seqlevels == ["chr1", "chr2", "chr3"]
    assert c.seqlengths == {"chr1": 100, "chr2": 200, "chr3": 300}


def test_concat_conflicting_seqinfo_raises() -> None:
    a = GRanges("chr1", IRanges(start=1, width=1), seqinfo=Seqinfo(["chr1"], [100]))
    b = GRanges("chr1", IRanges(start=1, width=1), seqinfo=Seqinfo(["chr1"], [200]))
    with pytest.raises(ValueError):
        a.concat(b)


def test_concat_remaps_codes_when_level_order_differs() -> None:
    a = GRanges(["x", "y"], IRanges(start=[1, 2], width=1))
    b = GRanges(["y", "x"], IRanges(start=[3, 4], width=1))
    c = a.concat(b)
    assert c.seqnames == ["x", "y", "y", "x"] and c.seqlevels == ["x", "y"]


def test_concat_mcols_must_match() -> None:
    a = GRanges("c", IRanges(start=1, width=1), mcols={"v": [1]})
    b = GRanges("c", IRanges(start=1, width=1))
    with pytest.raises(ValueError):
        a.concat(b)


def test_with_seqinfo_reorders_levels() -> None:
    x = gr().with_seqinfo(Seqinfo(["chr1", "chr2"], [50, 50]))
    assert x.seqnames == ["chr2", "chr1", "chr2"] and x.seqlevels == ["chr1", "chr2"]
    with pytest.raises(ValueError):
        gr().with_seqinfo(Seqinfo(["chr1"]))


def test_out_of_bound_warns() -> None:
    si = Seqinfo(["chr1", "chrC", "chrU"], [100, 100, None], [False, True, None])
    with pytest.warns(UserWarning, match="out-of-bound"):
        GRanges("chr1", IRanges(start=[90], end=[110]), seqinfo=si)
    with pytest.warns(UserWarning, match="out-of-bound"):
        GRanges("chr1", IRanges(start=[0], end=[10]), seqinfo=si)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        GRanges(["chrC", "chrU", "chr1"], IRanges(start=[90, 5000, 1], end=[110, 6000, 100]), seqinfo=si)


def test_trim() -> None:
    si = Seqinfo(["chr1", "chrC", "chrU"], [100, 100, None], [False, True, None])
    r = IRanges(start=[90, -5, 150, -20, 95], end=[110, 10, 160, -10, 99])
    with pytest.warns(UserWarning):
        x = GRanges(["chr1", "chr1", "chr1", "chr1", "chr1"], r, seqinfo=si)
    t = x.trim()
    assert list(zip(t.start.tolist(), t.end.tolist())) == [(90, 100), (1, 10), (101, 100), (1, 0), (95, 99)]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        t.trim()
    # circular and unknown-length sequences are left alone
    y = GRanges(["chrC", "chrU"], IRanges(start=[-5, 200], end=[500, 300]), seqinfo=si)
    assert y.trim().equals(y)


def test_equals_and_repr() -> None:
    assert gr().equals(gr())
    assert not gr().equals(gr().with_names(None))
    assert "GRanges of length 3" in repr(gr())
    assert "chr2:1-3 [+]" in repr(gr())


def test_empty() -> None:
    e = GRanges([], IRanges())
    assert len(e) == 0 and e.seqnames == [] and e.strand == []
    assert len(GRanges.combine([])) == 0
