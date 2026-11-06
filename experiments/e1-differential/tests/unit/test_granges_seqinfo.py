import warnings

import numpy as np
import polars as pl
import pytest

from ranges_e1 import GRanges, Hits, IRanges, Rle, Seqinfo


def test_seqlevels_follow_appearance_unless_seqinfo_given() -> None:
    g = GRanges(["chr2", "chr1", "chr2"], start=[1, 2, 3], width=1)
    assert g.seqlevels == ("chr2", "chr1")
    h = GRanges(["chr2"], start=[1], width=1, seqinfo=Seqinfo(["chr1", "chr2"]))
    assert h.seqlevels == ("chr1", "chr2")


def test_strand_validation_and_na_warning() -> None:
    with pytest.raises(ValueError):
        GRanges("chr1", start=1, width=1, strand="x")
    with pytest.warns(UserWarning):
        g = GRanges("chr1", start=1, width=1, strand=[None])
    assert g.strand.tolist() == ["*"]


def test_out_of_bound_warns_but_circular_and_unknown_do_not() -> None:
    with pytest.warns(UserWarning):
        GRanges("c", start=90, width=31, seqinfo=Seqinfo(["c"], [100]))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        GRanges("c", start=90, width=31, seqinfo=Seqinfo(["c"], [100], [True]))
        GRanges("c", start=90, width=31, seqinfo=Seqinfo(["c"]))


def test_seqinfo_merge_rules() -> None:
    a, b = Seqinfo(["x", "y"], [1, None]), Seqinfo(["y", "z"], [5, 7])
    with pytest.warns(UserWarning):
        m = a.merge(b)
    assert m.seqnames == ("x", "y", "z") and m.seqlengths == (1, 5, 7)
    with pytest.raises(ValueError):
        Seqinfo(["x"], [1]).merge(Seqinfo(["x"], [2]))
    with pytest.raises(ValueError):
        Seqinfo(["x"], genome="a").merge(Seqinfo(["x"], genome="b"))


def test_seqinfo_vector_lengths_are_strict_except_genome() -> None:
    assert Seqinfo(["a", "b"], genome="x").genome == ("x", "x")
    with pytest.raises(ValueError):
        Seqinfo(["a", "b"], seqlengths=5)


def test_trim_clips_and_leaves_unknown_lengths() -> None:
    with pytest.warns(UserWarning):
        g = GRanges(
            ["c", "c", "c", "u"], start=[90, 200, -5, -5], width=[31, 10, 3, 3],
            seqinfo=Seqinfo(["c", "u"], [100, None]),
        )
    t = g.trim()
    assert t.start.tolist() == [90, 101, 1, -5]
    assert t.width.tolist() == [11, 0, 0, 3]


def test_concat_merges_seqinfo_in_order() -> None:
    a = GRanges("c1", start=1, width=1, seqinfo=Seqinfo(["c1", "c2"], [10, 20]))
    b = GRanges("c3", start=1, width=1, seqinfo=Seqinfo(["c3", "c1"], [30, 10]))
    with pytest.warns(UserWarning):
        c = GRanges.concat([a, b])
    assert c.seqlevels == ("c1", "c2", "c3")
    assert c.seqnames.tolist() == ["c1", "c3"]


def test_setters_return_copies() -> None:
    g = GRanges("c", start=1, width=1)
    g2 = g.set_seqlengths([50]).set_names(["n"]).set_mcols(pl.DataFrame({"a": [1]}))
    assert g.seqlengths == {"c": None} and g2.seqlengths == {"c": 50} and g.names is None and g2.names == ("n",)


def test_hits_selects_and_counts() -> None:
    h = Hits(np.array([0, 0, 2]), np.array([3, 1, 2]), 4, 5)
    assert h.count_by_query().tolist() == [2, 0, 1, 0]
    first, last = h.select("first"), h.select("last")
    assert first.tolist() == [1, None, 2, None] and last.tolist() == [3, None, 2, None]
    with pytest.raises(ValueError):
        Hits(np.array([5]), np.array([0]), 2, 2)


def test_rle_roundtrip_merges_runs() -> None:
    r = Rle.encode(np.array([1, 1, 2, 2, 2, 1]))
    assert r.values.tolist() == [1, 2, 1] and r.lengths.tolist() == [2, 3, 1]
    assert r.decode().tolist() == [1, 1, 2, 2, 2, 1] and len(r) == 6


def test_na_results_are_masked_not_sentinels() -> None:
    q, s = IRanges(start=[1, 50], width=2), IRanges(start=[2], width=1)
    sel = q.find_overlaps(s, select="first")
    assert np.ma.isMaskedArray(sel) and sel.tolist() == [0, None]
    assert np.ma.isMaskedArray(q.nearest(IRanges(start=[], width=[])))


def test_recycled_coordinates_recycle_their_names_and_subset_works() -> None:
    g = GRanges(["a", "b"], IRanges([1], width=[2], names=["only"]))
    assert g.names == ("only", "only") and len(g) == 2
    assert g[1].names == ("only",)  # review finding 5: used to raise IndexError


def test_short_mcols_recycle_and_zero_rows_are_rejected() -> None:
    x = IRanges(start=[1, 3, 5], width=1)
    assert x.set_mcols(pl.DataFrame({"m": [7]})).mcols["m"].to_list() == [7, 7, 7]
    with pytest.warns(UserWarning):
        assert x.set_mcols(pl.DataFrame({"m": [7, 8]})).mcols["m"].to_list() == [7, 8, 7]
    with pytest.raises(ValueError):
        x.set_mcols(pl.DataFrame({"m": []}, schema={"m": pl.Int64}))
    with pytest.raises(ValueError):
        x.set_mcols(pl.DataFrame({"m": [1, 2, 3, 4]}))
