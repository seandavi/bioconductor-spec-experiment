"""Seqinfo (spec/90-seqinfo.md)."""

from __future__ import annotations

import pytest

from helpers import gr
from ranges_e3 import GRanges, IRanges, Seqinfo, concat


@pytest.mark.spec("SEQI.1")
def test_seqinfo_structure() -> None:
    si = Seqinfo(["b", "a", "c"], [100, 200, None], [False, True, None], ["g1", "g1", None])
    assert si.seqnames == ["b", "a", "c"]
    assert si.seqlengths == [100, 200, None]
    assert si.is_circular == [False, True, None]
    assert si.genome == ["g1", "g1", None]
    assert Seqinfo(["x"]).seqlengths == [None]
    with pytest.raises(ValueError):
        Seqinfo(["a", "a"])
    with pytest.raises(ValueError):
        Seqinfo(["a", "b"], seqlengths=[1])
    with pytest.raises(ValueError):
        Seqinfo(["a"], seqlengths=[-1])


@pytest.mark.spec("SEQI.2")
def test_merge_union_order_and_attribute_fill() -> None:
    x = Seqinfo(["a", "b"], [10, 20])
    y = Seqinfo(["c", "a"], [30, None])
    m = x.merge(y)
    assert m.seqnames == ["a", "b", "c"] and m.seqlengths == [10, 20, 30]
    circ = Seqinfo(["a"], is_circular=[True]).merge(Seqinfo(["a"]))
    assert circ.is_circular == [True]
    gen = Seqinfo(["a"]).merge(Seqinfo(["a"], genome=["hg19"]))
    assert gen.genome == ["hg19"]
    assert Seqinfo(["a"], [10]).merge(Seqinfo(["b"])).seqnames == ["a", "b"]
    assert Seqinfo(["a", "b"]).merge(Seqinfo(["b", "a"])).seqnames == ["a", "b"]


@pytest.mark.spec("SEQI.2")
@pytest.mark.parametrize(
    "x,y",
    [
        (Seqinfo(["a"], [10]), Seqinfo(["a"], [11])),
        (Seqinfo(["a"], is_circular=[True]), Seqinfo(["a"], is_circular=[False])),
        (Seqinfo(["a"], genome=["hg19"]), Seqinfo(["a"], genome=["hg38"])),
    ],
)
def test_merge_conflicts_refused(x: Seqinfo, y: Seqinfo) -> None:
    with pytest.raises(ValueError):
        x.merge(y)
    a = GRanges("a", IRanges(start=1, end=2), seqinfo=x)
    b = GRanges("a", IRanges(start=1, end=2), seqinfo=y)
    with pytest.raises(ValueError):
        concat(a, b)


@pytest.mark.spec("SEQI.3")
def test_granges_seqinfo_accessors() -> None:
    si = Seqinfo(["b", "a", "c"], [100, 200, None], [False, True, None], ["g1", "g1", None])
    g = GRanges(["b", "a", "b"], IRanges(start=[1, 2, 3], end=[3, 4, 5]), seqinfo=si)
    assert g.seqlevels == ["b", "a", "c"]
    assert g.seqlengths == [100, 200, None]
    assert g.is_circular == [False, True, None]
    assert g.genome == ["g1", "g1", None]
    assert g.seqinfo == si


@pytest.mark.spec("SEQI.4")
def test_default_seqinfo_and_unknown_seqname() -> None:
    g = gr(("chr2", 1, 2, "+"), ("chr1", 1, 2, "-"), ("chr10", 1, 2, "*"), ("chr1", 2, 3, "+"))
    assert g.seqlevels == ["chr2", "chr1", "chr10"]
    assert g.seqlengths == [None] * 3 and g.genome == [None] * 3
    with pytest.raises(ValueError):
        GRanges("z", IRanges(start=1, end=2), seqinfo=Seqinfo(["a"]))
    explicit = GRanges("a", IRanges(start=1, end=2), seqinfo=Seqinfo(["b", "a"]))
    assert explicit.seqlevels == ["b", "a"]


@pytest.mark.spec("SEQI.5")
def test_levels_survive_subsetting_and_ops() -> None:
    g = GRanges(["a", "b"], IRanges(start=[1, 2], end=[3, 4]), seqinfo=Seqinfo(["a", "b", "c"], [10, 20, 30]))
    assert g[0].seqlevels == ["a", "b", "c"] and g[0].seqlengths == [10, 20, 30]
    assert g[[]].seqlevels == ["a", "b", "c"]
    assert g.trim().seqinfo == g.seqinfo
