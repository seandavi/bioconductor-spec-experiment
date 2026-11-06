"""GRanges overlaps, neighbors, ordering (spec/95-granges.md GR.6, GR.14-16, GR.18)."""

from __future__ import annotations

import pytest

from helpers import gr, gs, na
from ranges_e3 import GRanges, IRanges, Seqinfo, concat


def pairs(q: GRanges, s: GRanges, **kw: object) -> list[tuple[int, int]]:
    return q.find_overlaps(s, **kw).pairs()  # type: ignore[union-attr]


@pytest.mark.spec("GR.6")
def test_strand_compatibility_via_overlaps() -> None:
    plus, minus, star = gr(("a", 1, 5, "+")), gr(("a", 1, 5, "-")), gr(("a", 1, 5, "*"))
    assert len(plus.find_overlaps(plus)) == 1 and len(plus.find_overlaps(star)) == 1
    assert len(star.find_overlaps(minus)) == 1 and len(plus.find_overlaps(minus)) == 0
    assert len(plus.find_overlaps(minus, ignore_strand=True)) == 1
    assert len(plus.find_overlaps(gr(("b", 1, 5, "+")), ignore_strand=True)) == 0


Q = gr(("a", 5, 10, "+"), ("a", 5, 10, "-"), ("a", 5, 10, "*"), ("b", 5, 10, "+"), ("a", 5, 10, "*"))
S = gr(("a", 1, 7, "+"), ("a", 1, 7, "-"), ("a", 1, 7, "*"), ("b", 1, 7, "-"))


@pytest.mark.spec("GR.14")
def test_overlaps_require_comparable_pairs() -> None:
    assert pairs(Q, S) == [(0, 0), (0, 2), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (4, 0), (4, 1), (4, 2)]
    assert pairs(Q, S, ignore_strand=True) == [
        (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (3, 3), (4, 0), (4, 1), (4, 2),
    ]
    assert Q.count_overlaps(S).tolist() == [2, 2, 3, 0, 3]
    assert Q.count_overlaps(S, ignore_strand=True).tolist() == [3, 3, 3, 1, 3]
    assert Q.overlaps_any(S).tolist() == [True, True, True, False, True]
    assert gs(Q.subset_by_overlaps(S)) == [gs(Q)[i] for i in (0, 1, 2, 4)]
    assert gs(Q.subset_by_overlaps(S, invert=True)) == [gs(Q)[3]]
    assert na(Q.find_overlaps(S, select="first")) == [0, 1, 0, None, 0]
    assert na(Q.find_overlaps(S, select="last")) == [2, 2, 2, None, 2]


@pytest.mark.spec("GR.14")
def test_overlaps_types_and_max_gap_within_comparable_pairs() -> None:
    eq = gr(("a", 1, 7, "+"), ("a", 1, 7, "-"))
    assert pairs(eq, S, type="equal") == [(0, 0), (0, 2), (1, 1), (1, 2)]
    assert pairs(eq, S, type="equal", ignore_strand=True) == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    assert pairs(gr(("a", 1, 5, "+")), gr(("a", 7, 9, "+")), max_gap=1) == [(0, 0)]
    assert pairs(gr(("a", 1, 5, "+")), gr(("a", 7, 9, "-")), max_gap=1) == []
    assert pairs(gr(("a", 6, 8, "+")), gr(("a", 1, 10, "+")), type="within") == [(0, 0)]
    with pytest.raises(ValueError):
        Q.find_overlaps(S, max_gap=0, min_overlap=1)


@pytest.mark.spec("GR.14")
def test_overlaps_match_seqnames_by_name_across_seqinfo() -> None:
    q = GRanges(["b", "a"], IRanges(start=1, end=5), seqinfo=Seqinfo(["b", "a"]))
    s = GRanges(["a", "b"], IRanges(start=1, end=5), seqinfo=Seqinfo(["a", "b"]))
    assert pairs(q, s) == [(0, 1), (1, 0)]
    assert len(gr(("a", 1, 5, "+")).find_overlaps(gr(("z", 1, 5, "+")))) == 0
    assert len(gr().find_overlaps(S)) == 0 and len(Q.find_overlaps(gr())) == 0


NS = gr(("a", 1, 3, "+"), ("a", 10, 12, "+"), ("a", 20, 22, "-"), ("a", 30, 32, "*"), ("b", 5, 7, "+"), ("a", 40, 42, "-"))
NQ = gr(("a", 6, 7, "+"), ("a", 6, 7, "-"), ("a", 6, 7, "*"), ("a", 6, 7, "+"), ("a", 25, 26, "+"), ("b", 1, 2, "-"),
        ("c", 1, 2, "+"), ("a", 50, 51, "*"))


@pytest.mark.spec("GR.15")
def test_nearest_restricted_to_comparable_subjects() -> None:
    assert na(NQ.nearest(NS)) == [0, 2, 0, 0, 3, None, None, 5]
    assert NQ.nearest(NS, select="all").pairs() == [
        (0, 0), (0, 1), (1, 2), (2, 0), (2, 1), (3, 0), (3, 1), (4, 3), (7, 5),
    ]
    assert na(NQ.nearest(NS, ignore_strand=True)) == [0, 0, 0, 0, 2, 4, None, 5]
    assert na(NS.nearest()) == [1, 0, 3, 2, None, 3]


@pytest.mark.spec("GR.15")
def test_precede_follow_direction_and_mirroring() -> None:
    assert na(NQ.precede(NS)) == [1, None, 1, 1, 3, None, None, 5]
    assert na(NQ.follow(NS)) == [0, 2, 0, 0, 1, None, None, 3]
    assert na(NQ.precede(NS, ignore_strand=True)) == [1, 1, 1, 1, 3, 4, None, None]
    assert na(NQ.follow(NS, ignore_strand=True)) == [0, 0, 0, 0, 2, None, None, 5]
    assert NQ.precede(NS, select="all").pairs() == [(0, 1), (2, 1), (3, 1), (4, 3), (7, 5)]
    assert NQ.follow(NS, select="all").pairs() == [(0, 0), (1, 2), (2, 0), (3, 0), (4, 1), (7, 3)]


@pytest.mark.spec("GR.15")
def test_star_query_merges_directions_by_distance() -> None:
    q = gr(("a", 100, 101, "*"))
    far_right_near_left = gr(("a", 140, 141, "+"), ("a", 90, 91, "-"))
    assert na(q.precede(far_right_near_left)) == [1]
    near_right_far_left = gr(("a", 110, 111, "+"), ("a", 70, 71, "-"))
    assert na(q.precede(near_right_far_left)) == [0]
    tie = gr(("a", 120, 121, "+"), ("a", 80, 81, "-"))
    assert q.precede(tie, select="all").pairs() == [(0, 0), (0, 1)]
    assert na(q.precede(tie)) == [0]
    only_minus_left = gr(("a", 80, 81, "-"), ("a", 70, 71, "-"))
    assert na(q.precede(only_minus_left)) == [0] and na(q.follow(only_minus_left)) == [None]
    minus_right = gr(("a", 120, 121, "-"), ("a", 130, 131, "-"))
    assert na(q.follow(minus_right)) == [0] and na(q.precede(minus_right)) == [None]


@pytest.mark.spec("GR.15", "NEAR.4")
def test_tie_rules_stay_with_the_operation_on_minus_queries() -> None:
    q = gr(("a", 6, 7, "-"))
    s = gr(("a", 1, 3, "-"), ("a", 1, 3, "-"), ("a", 10, 12, "-"), ("a", 10, 12, "-"))
    assert na(q.precede(s)) == [0] and na(q.follow(s)) == [3]


@pytest.mark.spec("GR.15")
def test_distance_is_na_for_non_comparable_pairs() -> None:
    q = gr(("a", 6, 7, "+"), ("a", 6, 7, "-"), ("a", 6, 7, "*"), ("a", 6, 7, "+"), ("a", 25, 26, "+"), ("b", 1, 2, "-"))
    s = gr(*[("a", 1, 3, "+")] * 6)
    assert na(q.distance(s)) == [2, None, 2, 2, 21, None]
    assert na(q.distance(s, ignore_strand=True)) == [2, 2, 2, 2, 21, None]
    d = gr(("a", 1, 3, "+"), ("a", 1, 3, "+"), ("a", 1, 3, "*")).distance(gr(("a", 5, 7, "-")))
    assert na(d) == [None, None, 1]
    with pytest.raises(ValueError):
        q.distance(s[[0, 1]])
    assert len(gr().distance(gr())) == 0


@pytest.mark.spec("GR.16")
def test_sort_order_is_level_then_strand_then_position() -> None:
    g = gr(("b", 5, 6, "+"), ("a", 1, 2, "-"), ("b", 3, 4, "-"), ("c", 1, 2, "*"), ("a", 1, 2, "+"), ("b", 3, 3, "+"),
           ("a", 1, 2, "+"), seqinfo=Seqinfo(["b", "a", "c"]))
    assert g.order().tolist() == [5, 0, 2, 4, 6, 1, 3]
    assert gs(g.sort()) == [
        ("b", 3, 3, "+"), ("b", 5, 6, "+"), ("b", 3, 4, "-"), ("a", 1, 2, "+"), ("a", 1, 2, "+"), ("a", 1, 2, "-"),
        ("c", 1, 2, "*"),
    ]
    assert g.order(decreasing=True).tolist() == [3, 1, 4, 6, 2, 0, 5]
    assert g.rank("first").tolist() == [2, 6, 3, 7, 4, 1, 5]
    assert g.duplicated().tolist() == [False] * 6 + [True]
    assert gs(g.unique()) == gs(g)[:6]


@pytest.mark.spec("GR.16")
def test_strand_beats_start_and_levels_are_not_alphabetical() -> None:
    t = gr(("a", 1, 2, "*"), ("a", 1, 2, "-"), ("a", 1, 2, "+"), ("a", 1, 2, "+"))
    assert t.order().tolist() == [2, 3, 1, 0]
    assert gr(("a", 1, 3, "+"), ("a", 2, 2, "+"), ("a", 1, 2, "+")).order().tolist() == [2, 0, 1]
    named = gr(("a", 3, 4, "+"), ("a", 1, 2, "+")).with_names(["x", "y"])
    assert named.sort().names == ["y", "x"]
    dup = concat(gr(("a", 1, 2, "+")), gr(("a", 1, 2, "-")), gr(("a", 1, 2, "+")))
    assert dup.duplicated().tolist() == [False, False, True] and dup.duplicated(from_last=True).tolist() == [True, False, False]


@pytest.mark.spec("GR.18")
def test_circular_levels_are_linear_but_not_trimmed() -> None:
    si = Seqinfo(["c", "n"], [30, 30], [True, False])
    g = GRanges(["c", "n"], IRanges(start=[25, 25], end=[35, 35]), seqinfo=si)
    assert [(a, b) for _, a, b, _ in gs(g.trim())] == [(25, 35), (25, 30)]
    assert g.is_circular == [True, False]
    other = GRanges("c", IRanges(start=1, end=4), seqinfo=Seqinfo(["c"], [30], [True]))
    assert len(g[[0]].find_overlaps(other)) == 0
    assert [len(r) for r in g.coverage().values()] == [30, 30]
    assert g.coverage()["c"].to_numpy()[:5].tolist() == [0, 0, 0, 0, 0]
