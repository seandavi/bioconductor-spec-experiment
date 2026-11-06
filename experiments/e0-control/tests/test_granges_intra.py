import pytest

from ranges_e0 import GRanges, IRanges, OutOfBoundWarning, Seqinfo

SI = Seqinfo(["chr1", "chr2"], [100, 50], [False, False], "g1")


def make() -> GRanges:
    return GRanges(
        ["chr1", "chr1", "chr1", "chr2"],
        IRanges([10, 10, 10, 5], width=[11, 11, 11, 2]),
        ["+", "-", "*", "-"],
        seqinfo=SI,
    )


def ivs(g):
    return [(s, e) for _, s, e, _ in g.to_tuples()]


def test_resize_minus_strand_anchors_at_end():
    g = make()
    # R: resize(g,4) -> 10-13 (+), 17-20 (-), 10-13 (*), chr2 5-6 (-) anchored at its end 6 -> 3-6
    assert ivs(g.resize(4)) == [(10, 13), (17, 20), (10, 13), (3, 6)]
    # fix="end" is the mirror image
    assert ivs(g.resize(4, fix="end")) == [(17, 20), (10, 13), (17, 20), (5, 8)]
    # R: center is strand independent
    assert ivs(g.resize(5, fix="center"))[:3] == [(13, 17)] * 3


def test_resize_ignore_strand():
    g = make()
    # R: resize(g, 4, ignore.strand=TRUE) treats everything as +
    assert ivs(g.resize(4, ignore_strand=True)) == [(10, 13), (10, 13), (10, 13), (5, 8)]


def test_flank_strand_aware():
    g = make()
    # R: flank(g,3): + upstream 7-9 ; - upstream is on the right 21-23 ; * like +
    assert ivs(g.flank(3))[:3] == [(7, 9), (21, 23), (7, 9)]
    assert ivs(g.flank(3, start=False))[:3] == [(21, 23), (7, 9), (21, 23)]
    assert ivs(g.flank(3, both=True))[:3] == [(7, 12), (18, 23), (7, 12)]
    assert ivs(g.flank(-3))[:3] == [(10, 12), (18, 20), (10, 12)]
    assert ivs(g.flank(3, ignore_strand=True))[:3] == [(7, 9)] * 3


def test_flank_zero_width_minus_strand():
    g = GRanges("chr2", IRanges(5, 4), "-", seqinfo=SI)
    assert ivs(g.flank(3)) == [(5, 7)]  # R: upstream of a - range is past its end
    assert ivs(g.flank(3, start=False)) == [(2, 4)]


def test_promoters_strand_aware():
    g = make()
    # R: promoters(g, 5, 3): + 5-12 ; - 18-25 ; * 5-12 ; chr2 (5-6 -): TSS = end 6 -> 4-11
    assert ivs(g.promoters(5, 3)) == [(5, 12), (18, 25), (5, 12), (4, 11)]


def test_shift_and_narrow_are_not_strand_aware():
    g = make()
    assert ivs(g.shift(2)) == [(12, 22), (12, 22), (12, 22), (7, 8)]
    # R: narrow(g, start=2, end=-2) -> 11-19 for every strand
    assert ivs(g.narrow(start=2, end=-2)) == [(11, 19), (11, 19), (11, 19), (6, 5)]
    assert ivs(g.narrow(start=2)) == [(11, 20), (11, 20), (11, 20), (6, 6)]
    assert ivs(g.narrow(end=-2)) == [(10, 19)] * 3 + [(5, 5)]


def test_intra_ops_keep_names_mcols_seqinfo_and_strand():
    g = make().with_names(list("abcd")).with_mcols({"v": [1, 2, 3, 4]})
    r = g.resize(3)
    assert r.names == list("abcd") and r.mcols["v"].to_list() == [1, 2, 3, 4]
    assert r.strand.tolist() == ["+", "-", "*", "-"] and r.seqinfo == SI
    assert g.resize(3, use_names=False).names is None


def test_shift_out_of_bounds_warns():
    g = make()
    with pytest.warns(OutOfBoundWarning):
        out = g.shift(95)
    assert ivs(out)[0] == (105, 115)  # R: warns, does not trim


def test_vector_arguments():
    g = make()
    assert ivs(g.resize([1, 2, 3, 1]))[:3] == [(10, 10), (19, 20), (10, 12)]
    assert ivs(g.shift([0, 1, 2, 3]))[0] == (10, 20)


def test_flank_start_argument_must_be_recyclable():
    g = make()
    with pytest.raises(ValueError, match="start"):
        g.flank(2, start=[True, False])
    assert len(g.flank(2, start=[True, False, True, False])) == 4
