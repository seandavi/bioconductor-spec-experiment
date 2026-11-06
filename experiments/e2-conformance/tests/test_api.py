"""Python-API behaviour the R-derived suite cannot see: 0-based indices, immutability, sentinels."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from ranges_e2 import NA_INDEX, GRanges, Hits, IRanges, Seqinfo


def test_subscripts_are_zero_based_numpy_style() -> None:
    x = IRanges(start=[1, 5, 12], width=[3, 4, 2], names=["a", "b", "c"])
    assert x[0].start.tolist() == [1]
    assert x[-1].start.tolist() == [12]
    assert x[[2, 0]].names == ["c", "a"]
    assert x[np.array([True, False, True])].start.tolist() == [1, 12]
    assert x[["b"]].start.tolist() == [5]
    assert x[1:].start.tolist() == [5, 12]
    with pytest.raises(IndexError):
        x[3]
    with pytest.raises(KeyError):
        x[["zzz"]]


def test_objects_are_immutable() -> None:
    x = IRanges(start=[1], width=[2])
    with pytest.raises(AttributeError):
        x.foo = 1  # type: ignore[attr-defined]
    g = GRanges(["a"], x)
    with pytest.raises(AttributeError):
        g.foo = 1  # type: ignore[attr-defined]


def test_ops_return_new_objects_and_leave_input_alone() -> None:
    x = IRanges(start=[10, 1], width=[3, 3], names=["a", "b"])
    y = x.shift(5)
    assert x.start.tolist() == [10, 1]
    assert y.start.tolist() == [15, 6]
    assert x.sort().names == ["b", "a"]
    assert x.names == ["a", "b"]


def test_select_uses_zero_based_indices_and_na_sentinel() -> None:
    q = IRanges(start=[1, 50], width=[5, 5])
    s = IRanges(start=[3, 100, 4], width=[5, 5, 5])
    assert q.find_overlaps(s, select="first").tolist() == [0, NA_INDEX]
    assert q.find_overlaps(s, select="last").tolist() == [2, NA_INDEX]
    hits = q.find_overlaps(s)
    assert isinstance(hits, Hits)
    assert hits.pairs() == [(0, 0), (0, 2)]
    assert (hits.query_length, hits.subject_length) == (2, 3)


def test_revmap_is_zero_based_list_column() -> None:
    x = IRanges(start=[1, 3, 20], width=[5, 5, 2])
    r = x.reduce(with_revmap=True)
    assert r.mcols["revmap"].to_list() == [[0, 1], [2]]


def test_mcols_without_columns_is_an_empty_frame() -> None:
    x = IRanges(start=[1, 2], width=[1, 1])
    assert not x.has_mcols
    assert x.mcols.width == 0
    y = x.with_mcols(pl.DataFrame({"score": [1.5, 2.5]}))
    assert y.has_mcols and y[[1]].mcols["score"].to_list() == [2.5]
    assert y.with_mcols(None).has_mcols is False


def test_int32_overflow_is_an_error() -> None:
    with pytest.raises(ValueError):
        IRanges(start=[2**31], width=[1])
    with pytest.raises(ValueError):
        IRanges(start=[2**31 - 2], width=[1]).shift(10)


def test_granges_accessors_and_seqlevel_order() -> None:
    si = Seqinfo(["chr2", "chr1"], [200, 100])
    g = GRanges(["chr1", "chr2", "chr1"], IRanges(start=[5, 1, 2], width=[1, 1, 1]), ["+", "-", "*"], seqinfo=si)
    assert g.seqlevels == ["chr2", "chr1"]
    assert g.sort().seqname_strings() == ["chr2", "chr1", "chr1"]
    assert g.seqnames.values.to_list() == ["chr1", "chr2", "chr1"]  # type: ignore[union-attr]
    assert g.strand.values.to_list() == ["+", "-", "*"]  # type: ignore[union-attr]
    assert g.seqlengths == (200, 100)


def test_granges_rejects_bad_input() -> None:
    r = IRanges(start=[1], width=[1])
    with pytest.raises(ValueError):
        GRanges(["a"], r, "x")
    with pytest.raises(ValueError):
        GRanges(["z"], r, seqinfo=Seqinfo(["a"]))
    with pytest.raises(ValueError):
        GRanges(["a", "b"], r)


def test_seqinfo_merge_conflict() -> None:
    with pytest.raises(ValueError):
        Seqinfo(["a"], [10]).merge(Seqinfo(["a"], [11]))
    merged = Seqinfo(["a"], [None], genome="hg").merge(Seqinfo(["a", "b"], [5, 6]))
    assert merged.seqnames == ("a", "b") and merged.seqlengths == (5, 6) and merged.genome == ("hg", None)
