"""GRanges construction, accessors, subsetting, concat, trim (spec/95-granges.md GR.1-6, GR.9)."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from helpers import gr, gs
from ranges_e3 import GRanges, IRanges, Seqinfo, concat


@pytest.mark.spec("GR.1")
def test_construction_recycling_and_defaults() -> None:
    g = GRanges(["a", "b"], IRanges(start=[1, 2], end=[3, 4]))
    assert gs(g) == [("a", 1, 3, "*"), ("b", 2, 4, "*")]
    g = GRanges("a", IRanges(start=[1, 2], end=[3, 4]), strand="-")
    assert gs(g) == [("a", 1, 3, "-"), ("a", 2, 4, "-")]
    g = GRanges(["a", "b"], IRanges(start=5, end=6), strand=["+", "-"])
    assert gs(g) == [("a", 5, 6, "+"), ("b", 5, 6, "-")]
    assert len(GRanges([], IRanges(start=[], width=[]))) == 0


@pytest.mark.spec("GR.1")
def test_construction_names_and_mcols() -> None:
    g = GRanges(["a", "b"], IRanges(start=[1, 2], end=[3, 4], names=["p", "q"]),
                mcols=pl.DataFrame({"score": [1.5, 2.0]}))
    assert g.names == ["p", "q"] and g.mcols["score"].to_list() == [1.5, 2.0]
    g2 = GRanges(["a", "b"], IRanges(start=[1, 2], end=[3, 4], names=["p", "q"]), names=["x", "y"])
    assert g2.names == ["x", "y"]
    assert GRanges("a", IRanges(start=1, end=2)).names is None


@pytest.mark.spec("GR.1")
def test_construction_length_mismatch_refused() -> None:
    with pytest.raises(ValueError):
        GRanges(["a", "b", "c"], IRanges(start=[1, 2], end=[3, 4]))
    with pytest.raises(ValueError):
        GRanges("a", IRanges(start=[1, 2], end=[3, 4]), strand=["+", "-", "*"])
    with pytest.raises(ValueError):
        GRanges("a", IRanges(start=[1, 2], end=[3, 4]), mcols=pl.DataFrame({"k": [1]}))


@pytest.mark.spec("GR.2")
def test_bad_strand_refused() -> None:
    with pytest.raises(ValueError):
        GRanges("a", IRanges(start=1, end=2), strand="x")
    with pytest.raises(ValueError):
        GRanges("a", IRanges(start=1, end=2), strand=["+", "."])


@pytest.mark.spec("GR.3")
def test_accessors() -> None:
    g = GRanges(["a", "b"], IRanges(start=[1, 2], end=[3, 4], names=["p", "q"]), strand=["+", "-"])
    assert g.seqnames.tolist() == ["a", "b"] and g.strand.tolist() == ["+", "-"]
    assert g.start.tolist() == [1, 2] and g.end.tolist() == [3, 4] and g.width.tolist() == [3, 3]
    assert isinstance(g.ranges, IRanges) and g.ranges.names is None and g.ranges.mcols.width == 0
    assert len(g) == 2 and g.mcols.width == 0


@pytest.mark.spec("GR.3")
def test_out_of_bounds_allowed_at_construction() -> None:
    g = GRanges("a", IRanges(start=[-5, 12], end=[3, 15]), seqinfo=Seqinfo(["a"], [10]))
    assert g.start.tolist() == [-5, 12]


@pytest.mark.spec("GR.4")
def test_subsetting_keeps_everything_but_selection() -> None:
    g = GRanges(["a", "b", "a"], IRanges(start=[1, 2, 3], end=[3, 4, 5], names=["p", "q", "r"]),
                strand=["+", "-", "*"], mcols=pl.DataFrame({"k": [1, 2, 3]}))
    s = g[[2, 0]]
    assert gs(s) == [("a", 3, 5, "*"), ("a", 1, 3, "+")]
    assert s.names == ["r", "p"] and s.mcols["k"].to_list() == [3, 1]
    assert gs(g["q"]) == [("b", 2, 4, "-")]
    assert gs(g[np.array([True, False, True])]) == [("a", 1, 3, "+"), ("a", 3, 5, "*")]
    assert gs(g[-1]) == [("a", 3, 5, "*")]
    assert g[[1]].seqlevels == ["a", "b"]
    with pytest.raises(IndexError):
        g[5]


@pytest.mark.spec("GR.5")
def test_concat_merges_seqinfo_and_remaps_codes() -> None:
    a = gr(("chr2", 1, 2, "+"))
    b = gr(("chr3", 1, 2, "-"))
    c = concat(a, b)
    assert gs(c) == [("chr2", 1, 2, "+"), ("chr3", 1, 2, "-")] and c.seqlevels == ["chr2", "chr3"]
    assert concat(b, a).seqlevels == ["chr3", "chr2"]
    c1 = GRanges("a", IRanges(start=1, end=2), seqinfo=Seqinfo(["a", "b"], [10, 20]))
    c2 = GRanges("c", IRanges(start=1, end=2), seqinfo=Seqinfo(["c", "a"], [30, None]))
    cc = concat(c1, c2)
    assert cc.seqlevels == ["a", "b", "c"] and cc.seqlengths == [10, 20, 30]
    assert gs(cc) == [("a", 1, 2, "*"), ("c", 1, 2, "*")]
    swapped = concat(c2, c1)
    assert swapped.seqlevels == ["c", "a", "b"] and gs(swapped) == [("c", 1, 2, "*"), ("a", 1, 2, "*")]


@pytest.mark.spec("GR.5")
def test_concat_names_and_mcols() -> None:
    a = GRanges("a", IRanges(start=1, end=2), names=["p"], mcols=pl.DataFrame({"x": [1]}))
    b = GRanges("a", IRanges(start=3, end=4), mcols=pl.DataFrame({"y": ["q"]}))
    c = concat(a, b)
    assert c.names == ["p", ""]
    assert c.mcols["x"].to_list() == [1, None] and c.mcols["y"].to_list() == [None, "q"]
    with pytest.raises(TypeError):
        concat(a, IRanges(start=1, end=2))  # type: ignore[type-var]


@pytest.mark.spec("GR.9")
def test_trim() -> None:
    si = Seqinfo(["a", "b", "c"], [10, None, 10], [False, None, True])
    g = GRanges(["a", "a", "a", "a", "b", "c"],
                IRanges(start=[12, -5, -5, 1, -3, -3], end=[15, -3, 200, 0, 99, 99]), seqinfo=si)
    t = g.trim()
    assert [(a, b) for _, a, b, _ in gs(t)] == [(11, 10), (1, 0), (1, 10), (1, 0), (-3, 99), (-3, 99)]
    z = GRanges("a", IRanges(start=[11, 12], end=[10, 11]), seqinfo=Seqinfo(["a"], [10])).trim()
    assert [(a, b) for _, a, b, _ in gs(z)] == [(11, 10), (11, 10)]
    inside = gr(("a", 2, 5, "+"), seqinfo=Seqinfo(["a"], [10]), names=["n"])
    assert gs(inside.trim()) == [("a", 2, 5, "+")] and inside.trim().names == ["n"]
