"""Strand-aware intra-range operations on GRanges (spec/95-granges.md GR.7-8)."""

from __future__ import annotations

import polars as pl
import pytest

from helpers import gr, gs
from ranges_e3 import GRanges, IRanges, Seqinfo

G = gr(("a", 10, 14, "+"), ("a", 10, 14, "-"), ("a", 10, 14, "*"), ("b", 5, 7, "-"),
       seqinfo=Seqinfo(["a", "b"], [100, 50]))


def span(x: GRanges) -> list[tuple[int, int]]:
    return [(s, e) for _, s, e, _ in gs(x)]


@pytest.mark.spec("GR.7")
def test_shift_and_narrow_ignore_strand() -> None:
    assert span(G.shift(2)) == [(12, 16)] * 3 + [(7, 9)]
    assert span(G.narrow(start=2, width=2)) == [(11, 12)] * 3 + [(6, 7)]
    assert span(G.narrow(end=3)) == [(10, 12)] * 3 + [(5, 7)]
    assert gs(G.shift(1))[1][3] == "-"


@pytest.mark.spec("GR.8")
def test_resize_minus_strand_swaps_fix() -> None:
    assert span(G.resize(2)) == [(10, 11), (13, 14), (10, 11), (6, 7)]
    assert span(G.resize(2, fix="end")) == [(13, 14), (10, 11), (13, 14), (5, 6)]
    assert span(G.resize(2, fix="center")) == [(11, 12)] * 3 + [(5, 6)]
    assert span(G.resize(8, fix="center")) == [(8, 15)] * 3 + [(2, 9)]
    assert span(G.resize(2, ignore_strand=True)) == [(10, 11)] * 3 + [(5, 6)]
    assert span(G.resize([1, 2, 3, 4])) == [(10, 10), (13, 14), (10, 12), (4, 7)]


@pytest.mark.spec("GR.8")
def test_flank_minus_strand_mirrors() -> None:
    assert span(G.flank(3)) == [(7, 9), (15, 17), (7, 9), (8, 10)]
    assert span(G.flank(3, start=False)) == [(15, 17), (7, 9), (15, 17), (2, 4)]
    assert span(G.flank(-3)) == [(10, 12), (12, 14), (10, 12), (5, 7)]
    assert span(G.flank(3, both=True)) == [(7, 12), (12, 17), (7, 12), (5, 10)]
    assert span(G.flank(3, start=False, both=True)) == [(12, 17), (7, 12), (12, 17), (2, 7)]
    assert span(G.flank(3, ignore_strand=True)) == [(7, 9), (7, 9), (7, 9), (2, 4)]
    assert span(G.flank(3, start=[True, False, True, True])) == [(7, 9), (7, 9), (7, 9), (8, 10)]


@pytest.mark.spec("GR.8")
def test_promoters_minus_strand() -> None:
    assert span(G.promoters(4, 2)) == [(6, 11), (13, 18), (6, 11), (6, 11)]
    assert span(G.promoters()) == [(-1990, 209), (-185, 2014), (-1990, 209), (-192, 2007)]
    assert span(G.promoters(4, 2, ignore_strand=True)) == [(6, 11)] * 3 + [(1, 6)]
    with pytest.raises(ValueError):
        G.promoters(-1, 0)


@pytest.mark.spec("GR.7", "GR.8", "CONV.8")
def test_intra_keeps_seqinfo_names_mcols_strand() -> None:
    g = G.with_names(["p", "q", "r", "s"]).with_mcols(pl.DataFrame({"k": [1, 2, 3, 4]}))
    for out in (g.shift(1), g.narrow(start=1), g.resize(1), g.flank(1), g.promoters(1, 1)):
        assert out.names == ["p", "q", "r", "s"] and out.mcols["k"].to_list() == [1, 2, 3, 4]
        assert out.seqinfo == G.seqinfo and out.strand.tolist() == ["+", "-", "*", "-"]
        assert out.seqnames.tolist() == ["a", "a", "a", "b"]
    empty = GRanges([], IRanges(start=[], width=[]))
    assert len(empty.resize(3)) == 0 and len(empty.flank(3)) == 0
