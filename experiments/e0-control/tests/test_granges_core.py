import numpy as np
import polars as pl
import pytest

from ranges_e0 import GRanges, IRanges, OutOfBoundWarning, Seqinfo, concat


def gr(spec, **kw):
    """Build a GRanges from 'chr:start-end:strand' strings (test helper)."""
    seq, start, end, strand = [], [], [], []
    for item in spec:
        c, rest, s = item.split(":")
        a, b = rest.split("-")
        seq.append(c)
        start.append(int(a))
        end.append(int(b))
        strand.append(s)
    return GRanges(seq, IRanges(start, end), strand, **kw)


def test_construction_and_accessors():
    g = GRanges(["chr2", "chr10", "chr1", "chr2"], IRanges([5, 3, 9, 1], width=[2, 3, 0, 4]), ["-", "+", "*", "+"])
    # R: seqlevels follow first appearance, not sorting
    assert g.seqlevels == ("chr2", "chr10", "chr1")
    assert g.seqnames.tolist() == ["chr2", "chr10", "chr1", "chr2"]
    assert g.strand.tolist() == ["-", "+", "*", "+"]
    assert g.start.tolist() == [5, 3, 9, 1] and g.end.tolist() == [6, 5, 8, 4]
    assert g.width.tolist() == [2, 3, 0, 4]
    assert g.ranges.to_tuples() == [(5, 6), (3, 5), (9, 8), (1, 4)]
    assert g.names is None and len(g) == 4
    assert g.seqlengths == {"chr2": None, "chr10": None, "chr1": None}


def test_default_strand_is_star_and_scalars_recycle():
    g = GRanges("chr1", IRanges([1, 2], [3, 4]))
    assert g.strand.tolist() == ["*", "*"]
    assert GRanges(["a", "b"], IRanges(1, 2), "+").strand.tolist() == ["+", "+"]


def test_construction_errors():
    with pytest.raises(ValueError, match="must have the length of the object"):
        GRanges(["a", "b", "c"], IRanges([1, 2], [3, 4]))
    with pytest.raises(ValueError, match="strand values"):
        GRanges("a", IRanges(1, 2), "x")
    with pytest.raises(ValueError, match="no entries in 'seqinfo'"):
        GRanges("chr1", IRanges(1, 2), seqinfo=Seqinfo(["chr2"]))
    with pytest.raises(ValueError):
        GRanges("a", IRanges(1, 2), seqinfo=Seqinfo(["a"]), seqlengths={"a": 3})


def test_seqlengths_argument_defines_levels():
    # R: seqlengths=c(chr1=10, chr2=20, chr3=30) adds chr3 although unused
    g = GRanges(["chr1", "chr2"], IRanges(1, 2), seqlengths={"chr1": 10, "chr2": 20, "chr3": 30})
    assert g.seqlevels == ("chr1", "chr2", "chr3")
    assert g.seqlengths["chr3"] == 30


def test_names_come_from_ranges_or_argument():
    g = GRanges("chr1", IRanges([1, 2], [3, 4], names=["x", "y"]))
    assert g.names == ["x", "y"]
    assert GRanges("chr1", IRanges([1, 2], [3, 4]), names=["p", "q"]).names == ["p", "q"]


def test_mcols_and_subsetting():
    g = GRanges(["a", "b", "a"], IRanges([1, 5, 9], width=2), ["+", "-", "*"], mcols={"score": [1.0, 2.0, 3.0]},
                names=["r1", "r2", "r3"])
    s = g[[2, 0]]
    assert s.to_tuples() == [("a", 9, 10, "*"), ("a", 1, 2, "+")]
    assert s.mcols["score"].to_list() == [3.0, 1.0]
    assert s.names == ["r3", "r1"]
    assert g["r2"].to_tuples() == [("b", 5, 6, "-")]
    assert g[np.array([True, False, True])].seqnames.tolist() == ["a", "a"]
    assert g[g.strand == "-"].to_tuples() == [("b", 5, 6, "-")]
    assert g[-1].to_tuples() == [("a", 9, 10, "*")]
    # R: subsetting keeps the full seqinfo even when levels become unused
    assert g[[0]].seqlevels == ("a", "b")
    assert len(g[[]]) == 0 and g[[]].seqlevels == ("a", "b")


def test_with_helpers_and_equality():
    g = gr(["chr1:1-5:+", "chr1:3-4:-"])
    assert g.with_strand("*").strand.tolist() == ["*", "*"]
    assert g.with_names(["a", "b"]).names == ["a", "b"]
    assert g == gr(["chr1:1-5:+", "chr1:3-4:-"])
    assert g != g.with_strand("*")
    assert g != g.with_seqinfo(Seqinfo(["chr1", "chr2"]))


def test_with_seqinfo_recodes_by_name():
    g = gr(["chr1:1-5:+", "chr2:3-4:-"])
    h = g.with_seqinfo(Seqinfo(["chr2", "chrX", "chr1"], [10, None, 20]))
    assert h.seqnames.tolist() == ["chr1", "chr2"]
    assert h.seqlengths["chr1"] == 20
    with pytest.raises(ValueError, match="absent"):
        g.with_seqinfo(Seqinfo(["chr1"]))


def test_out_of_bound_warning_on_construction():
    si = Seqinfo(["chr1"], [100])
    # R warns (does not fail) for ranges past the sequence end
    with pytest.warns(OutOfBoundWarning, match="1 out-of-bound"):
        GRanges("chr1", IRanges(1, 150), seqinfo=si)
    with pytest.warns(OutOfBoundWarning):
        GRanges("chr1", IRanges(-3, 5), seqinfo=si)


def test_no_out_of_bound_warning_for_circular_or_unknown_length(recwarn):
    GRanges("c1", IRanges(90, 120), seqinfo=Seqinfo(["c1"], [100], [True]))
    GRanges("c1", IRanges(90, 120), seqinfo=Seqinfo(["c1"], [None]))
    GRanges("c1", IRanges(101, 100), seqinfo=Seqinfo(["c1"], [100]))  # empty range at len+1
    assert len(recwarn) == 0


def test_trim():
    si = Seqinfo(["chr1", "chr2"], [100, 50], [False, False])
    with pytest.warns(OutOfBoundWarning):
        g = GRanges("chr1", IRanges([-5, 90, 150, 101, -5], [3, 110, 160, 101, -2]), seqinfo=si)
    # R: trim -> 1-3, 90-100, (101,100), (101,100), (1,0)
    assert g.trim().ranges.to_tuples() == [(1, 3), (90, 100), (101, 100), (101, 100), (1, 0)]
    assert g.trim().start.tolist() == [1, 90, 101, 101, 1]


def test_trim_skips_circular_and_unknown_lengths():
    si = Seqinfo(["c", "u", "k"], [100, None, 10], [True, False, False])
    # only "k" (known length, non-circular) is clipped; the warning applies to k as well
    with pytest.warns(OutOfBoundWarning):
        g = GRanges(["c", "u", "k"], IRanges([90, 90, 5], [120, 120, 15]), seqinfo=si)
    assert g.trim().ranges.to_tuples() == [(90, 120), (90, 120), (5, 10)]


def test_concat_merges_seqinfo_names_and_mcols():
    a = GRanges("chr1", IRanges(1, 5), "+", seqinfo=Seqinfo(["chr1", "chr2"], [100, None]), mcols={"x": [1.0]})
    b = GRanges("chr3", IRanges(2, 6), "-", seqinfo=Seqinfo(["chr3", "chr2", "chr1"], [30, 50, None]), mcols={"y": ["u"]})
    c = concat(a, b)
    # R: seqlevels chr1 chr2 chr3 ; lengths 100 50 30 ; mcols aligned with NA
    assert c.seqlevels == ("chr1", "chr2", "chr3")
    assert c.seqlengths == {"chr1": 100, "chr2": 50, "chr3": 30}
    assert c.to_tuples() == [("chr1", 1, 5, "+"), ("chr3", 2, 6, "-")]
    assert c.mcols["x"].to_list() == [1.0, None] and c.mcols["y"].to_list() == [None, "u"]
    assert GRanges.concat(a, a).strand.tolist() == ["+", "+"]


def test_concat_seqinfo_conflicts():
    a = GRanges("chr1", IRanges(1, 5), seqinfo=Seqinfo(["chr1"], [100]))
    b = GRanges("chr1", IRanges(1, 5), seqinfo=Seqinfo(["chr1"], [7]))
    with pytest.raises(ValueError, match="incompatible seqlengths"):
        concat(a, b)
    # circularity NA merges with a value (R: result has 1 circular sequence)
    c = GRanges("chr1", IRanges(1, 5), seqinfo=Seqinfo(["chr1"], None, [True]))
    d = GRanges("chr1", IRanges(1, 5), seqinfo=Seqinfo(["chr1"], [100]))
    m = concat(d, c)
    assert m.is_circular == {"chr1": True} and m.seqlengths == {"chr1": 100}


def test_concat_keeps_strand_and_names():
    a = gr(["chr1:1-2:+"]).with_names(["a"])
    b = gr(["chr1:3-4:-"])
    c = concat(a, b)
    assert c.names == ["a", ""] and c.strand.tolist() == ["+", "-"]


def test_start_is_read_only():
    g = gr(["chr1:1-5:+"])
    with pytest.raises(ValueError):
        g.start[0] = 3


def test_empty():
    e = GRanges.empty()
    assert len(e) == 0 and e.seqnames.tolist() == [] and len(e.sort()) == 0
