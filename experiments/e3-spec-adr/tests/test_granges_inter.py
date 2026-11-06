"""GRanges grouping ops, set ops, parallel ops, coverage (spec/95-granges.md GR.10-13, GR.17)."""

from __future__ import annotations

import polars as pl
import pytest

from helpers import gr, gs
from ranges_e3 import GRanges, IRanges, Seqinfo

R = gr(("b", 1, 3, "+"), ("a", 5, 6, "+"), ("a", 7, 9, "+"), ("a", 5, 6, "-"), ("a", 6, 8, "*"), ("b", 2, 6, "+"),
       seqinfo=Seqinfo(["b", "a", "c"], [20, 30, 40]))


@pytest.mark.spec("GR.10")
def test_reduce_groups_by_seqname_and_strand() -> None:
    assert gs(R.reduce()) == [("b", 1, 6, "+"), ("a", 5, 9, "+"), ("a", 5, 6, "-"), ("a", 6, 8, "*")]
    assert gs(R.reduce(ignore_strand=True)) == [("b", 1, 6, "*"), ("a", 5, 9, "*")]
    assert gs(R.reduce(min_gapwidth=2)) == gs(R.reduce())
    assert R.reduce().seqinfo == R.seqinfo
    assert gs(gr(("a", 1, 3, "+"), ("a", 5, 4, "+")).reduce(drop_empty_ranges=True)) == [("a", 1, 3, "+")]


@pytest.mark.spec("GR.10")
def test_reduce_output_order_is_level_then_strand_then_start() -> None:
    g = gr(("a", 9, 10, "*"), ("a", 5, 6, "-"), ("a", 1, 2, "+"), ("b", 1, 2, "+"), seqinfo=Seqinfo(["b", "a"]))
    assert gs(g.reduce()) == [("b", 1, 2, "+"), ("a", 1, 2, "+"), ("a", 5, 6, "-"), ("a", 9, 10, "*")]


@pytest.mark.spec("GR.10")
def test_revmap_refers_to_input_indices() -> None:
    r = R.reduce(with_revmap=True)
    assert r.mcols["revmap"].to_list() == [[0, 5], [1, 2], [3], [4]]
    ri = R.reduce(with_revmap=True, ignore_strand=True)
    assert ri.mcols["revmap"].to_list() == [[0, 5], [1, 3, 4, 2]]
    d = gr(("a", 1, 5, "+"), ("a", 3, 8, "+"), ("a", 1, 5, "-")).disjoin(with_revmap=True)
    assert gs(d) == [("a", 1, 2, "+"), ("a", 3, 5, "+"), ("a", 6, 8, "+"), ("a", 1, 5, "-")]
    assert d.mcols["revmap"].to_list() == [[0], [0, 1], [1], [2]]
    assert R.reduce(drop_empty_ranges=True, with_revmap=True).mcols["revmap"].to_list() == [[0, 5], [1, 2], [3], [4]]


@pytest.mark.spec("GR.10")
def test_range_and_disjoin() -> None:
    assert gs(R.range()) == [("b", 1, 6, "+"), ("a", 5, 9, "+"), ("a", 5, 6, "-"), ("a", 6, 8, "*")]
    assert gs(R.range(ignore_strand=True)) == [("b", 1, 6, "*"), ("a", 5, 9, "*")]
    assert gs(R.disjoin()) == [
        ("b", 1, 1, "+"), ("b", 2, 3, "+"), ("b", 4, 6, "+"), ("a", 5, 6, "+"), ("a", 7, 9, "+"),
        ("a", 5, 6, "-"), ("a", 6, 8, "*"),
    ]
    assert gs(R.disjoin(ignore_strand=True)) == [
        ("b", 1, 1, "*"), ("b", 2, 3, "*"), ("b", 4, 6, "*"), ("a", 5, 5, "*"), ("a", 6, 6, "*"),
        ("a", 7, 8, "*"), ("a", 9, 9, "*"),
    ]
    assert len(gr().reduce()) == 0 and len(gr().range()) == 0 and len(gr().disjoin()) == 0


@pytest.mark.spec("GR.11")
def test_gaps_with_known_lengths_emit_all_strands() -> None:
    assert gs(R.gaps()) == [
        ("b", 7, 20, "+"), ("b", 1, 20, "-"), ("b", 1, 20, "*"),
        ("a", 1, 4, "+"), ("a", 10, 30, "+"), ("a", 1, 4, "-"), ("a", 7, 30, "-"), ("a", 1, 5, "*"), ("a", 9, 30, "*"),
        ("c", 1, 40, "+"), ("c", 1, 40, "-"), ("c", 1, 40, "*"),
    ]
    assert gs(R.gaps(ignore_strand=True)) == [("b", 7, 20, "*"), ("a", 1, 4, "*"), ("a", 10, 30, "*"), ("c", 1, 40, "*")]
    assert gs(R.gaps(start=1, end=10))[:3] == [("b", 7, 10, "+"), ("b", 1, 10, "-"), ("b", 1, 10, "*")]


@pytest.mark.spec("GR.11")
def test_gaps_with_unknown_lengths() -> None:
    g = gr(("a", 5, 10, "+"), ("a", 20, 25, "+"), ("b", 3, 4, "-"), ("c", 1, 2, "*"),
           seqinfo=Seqinfo(["a", "b", "c", "d"], [None, 30, None, None]))
    assert gs(g.gaps()) == [
        ("a", 1, 4, "+"), ("a", 11, 19, "+"), ("b", 1, 30, "+"), ("b", 1, 2, "-"), ("b", 5, 30, "-"), ("b", 1, 30, "*"),
    ]
    assert gs(g.gaps(ignore_strand=True)) == [("a", 1, 4, "*"), ("a", 11, 19, "*"), ("b", 1, 2, "*"), ("b", 5, 30, "*")]
    forced = g.gaps(start=1, end=40)
    assert ("d", 1, 40, "+") in gs(forced) and ("c", 3, 40, "*") in gs(forced) and ("a", 26, 40, "+") in gs(forced)


A = gr(("a", 1, 8, "+"), ("a", 10, 20, "-"), ("b", 1, 5, "+"), seqinfo=Seqinfo(["a", "b"]))
B = gr(("a", 5, 12, "+"), ("a", 10, 15, "*"), ("b", 3, 9, "-"), seqinfo=Seqinfo(["a", "b"]))


@pytest.mark.spec("GR.12")
def test_set_operations_treat_star_as_own_strand() -> None:
    assert gs(A.union(B)) == [("a", 1, 12, "+"), ("a", 10, 20, "-"), ("a", 10, 15, "*"), ("b", 1, 5, "+"), ("b", 3, 9, "-")]
    assert gs(A.intersect(B)) == [("a", 5, 8, "+")]
    assert gs(A.setdiff(B)) == [("a", 1, 4, "+"), ("a", 10, 20, "-"), ("b", 1, 5, "+")]
    assert gs(A.union(B, ignore_strand=True)) == [("a", 1, 20, "*"), ("b", 1, 9, "*")]
    assert gs(A.intersect(B, ignore_strand=True)) == [("a", 5, 8, "*"), ("a", 10, 15, "*"), ("b", 3, 5, "*")]
    assert gs(A.setdiff(B, ignore_strand=True)) == [("a", 1, 4, "*"), ("a", 16, 20, "*"), ("b", 1, 2, "*")]


@pytest.mark.spec("GR.12")
def test_set_operations_merge_seqinfo() -> None:
    x = gr(("a", 1, 5, "+"), seqinfo=Seqinfo(["a"], [10]))
    y = gr(("b", 1, 5, "+"), seqinfo=Seqinfo(["b", "a"], [None, 10]))
    u = x.union(y)
    assert u.seqlevels == ["a", "b"] and u.seqlengths == [10, None]
    assert gs(u) == [("a", 1, 5, "+"), ("b", 1, 5, "+")]
    assert gs(x.intersect(y)) == [] and gs(x.setdiff(y)) == [("a", 1, 5, "+")]


@pytest.mark.spec("GR.13")
def test_parallel_ops_result_strand_and_ranges() -> None:
    x = gr(("a", 1, 8, "+"))
    assert gs(x.pintersect(gr(("a", 5, 12, "+")))) == [("a", 5, 8, "+")]
    assert gs(x.pintersect(gr(("a", 5, 12, "*")))) == [("a", 5, 8, "+")]
    assert gs(gr(("a", 1, 8, "*")).pintersect(gr(("a", 5, 12, "-")))) == [("a", 5, 8, "-")]
    assert gs(x.punion(gr(("a", 5, 12, "*")))) == [("a", 1, 12, "+")]
    assert gs(x.psetdiff(gr(("a", 5, 12, "*")))) == [("a", 1, 4, "+")]
    assert gs(gr(("a", 1, 8, "*")).psetdiff(gr(("a", 5, 12, "+")))) == [("a", 1, 4, "+")]


@pytest.mark.spec("GR.13")
def test_parallel_ops_non_comparable_pairs() -> None:
    x = gr(("a", 1, 8, "+"))
    minus, other = gr(("a", 5, 12, "-")), gr(("b", 5, 12, "+"))
    assert gs(x.pintersect(minus)) == [("a", 1, 0, "+")]
    assert gs(x.pintersect(other)) == [("a", 1, 0, "+")]
    assert gs(x.psetdiff(minus)) == [("a", 1, 8, "+")]
    assert gs(x.psetdiff(other)) == [("a", 1, 8, "+")]
    with pytest.raises(ValueError):
        x.punion(minus)
    with pytest.raises(ValueError):
        x.punion(gr(("a", 1, 2, "+"), ("a", 1, 2, "+")))
    assert gs(x.pintersect(minus, ignore_strand=True)) == [("a", 5, 8, "+")]
    assert gs(gr(("a", 1, 8, "*")).pintersect(minus, ignore_strand=True)) == [("a", 5, 8, "*")]
    assert gs(gr(("a", 1, 8, "-")).punion(gr(("a", 5, 12, "+")), ignore_strand=True)) == [("a", 1, 12, "-")]
    assert gs(gr(("a", 1, 8, "-")).psetdiff(gr(("a", 5, 12, "+")), ignore_strand=True)) == [("a", 1, 4, "-")]
    with pytest.raises(ValueError):
        x.punion(other, ignore_strand=True)


@pytest.mark.spec("GR.13", "CONV.8")
def test_parallel_keeps_x_names_only() -> None:
    x = gr(("a", 1, 8, "+")).with_names(["n"]).with_mcols(pl.DataFrame({"k": [1]}))
    out = x.pintersect(gr(("a", 5, 12, "+")))
    assert out.names == ["n"] and out.mcols.width == 0


@pytest.mark.spec("GR.17")
def test_coverage_per_level_ignores_strand() -> None:
    g = gr(("a", 1, 4, "+"), ("a", 3, 6, "-"), ("b", 2, 3, "+"), ("a", 10, 12, "+"),
           seqinfo=Seqinfo(["a", "b", "c"], [15, 5, 4]))
    cov = g.coverage()
    assert list(cov) == ["a", "b", "c"]
    runs = lambda r: list(zip(r.lengths.tolist(), r.values.tolist()))  # noqa: E731
    assert runs(cov["a"]) == [(2, 1), (2, 2), (2, 1), (3, 0), (3, 1), (3, 0)]
    assert runs(cov["b"]) == [(1, 0), (2, 1), (2, 0)]
    assert runs(cov["c"]) == [(4, 0)]
    assert runs(g.coverage(shift=1)["a"]) == [(1, 0), (2, 1), (2, 2), (2, 1), (3, 0), (3, 1), (2, 0)]


@pytest.mark.spec("GR.17")
def test_coverage_length_rules() -> None:
    g = gr(("a", 1, 4, "+"), ("a", 3, 6, "+"))
    assert [(a, b) for a, b in zip(g.coverage()["a"].lengths, g.coverage()["a"].values)] == [(2, 1), (2, 2), (2, 1)]
    assert len(g.coverage(width=10)["a"]) == 10
    short = GRanges(["a", "b"], IRanges(start=1, end=3), seqinfo=Seqinfo(["a", "b"], [2, None]))
    cov = short.coverage()
    assert len(cov["a"]) == 2 and len(cov["b"]) == 3
    assert len(GRanges("a", IRanges(start=1, end=3), seqinfo=Seqinfo(["a", "z"])).coverage()["z"]) == 0
