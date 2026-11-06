import numpy as np
import pytest

from ranges_e0 import NO_MATCH, GRanges, IRanges, Seqinfo

NA = NO_MATCH


def pairs(h):
    return [tuple(p) for p in h.to_numpy().tolist()]


def Q():
    return GRanges(["chr1", "chr1", "chr1", "chr2", "chr1"], IRanges([1, 1, 1, 1, 50], width=10),
                   ["+", "-", "*", "+", "+"])


def S():
    return GRanges(["chr1", "chr1", "chr1", "chr2"], IRanges([5, 5, 5, 5], width=3), ["+", "-", "*", "-"])


def test_find_overlaps_strand_compatibility():
    # R (1-based): (1,1) (1,3) (2,2) (2,3) (3,1) (3,2) (3,3)  - '*' matches any strand
    assert pairs(Q().find_overlaps(S())) == [(0, 0), (0, 2), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]


def test_find_overlaps_ignore_strand():
    # R: adds (1,2) (4,4)... : every same-seqname pair
    assert pairs(Q().find_overlaps(S(), ignore_strand=True)) == [
        (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (3, 3)
    ]


def test_select_count_any_subset():
    # R: select="first" -> 1 2 1 NA NA ; countOverlaps -> 2 2 3 0 0
    assert Q().find_overlaps(S(), select="first").tolist() == [0, 1, 0, NA, NA]
    assert Q().find_overlaps(S(), select="last").tolist() == [2, 2, 2, NA, NA]
    assert Q().count_overlaps(S()).tolist() == [2, 2, 3, 0, 0]
    assert Q().overlaps_any(S()).tolist() == [True, True, True, False, False]
    assert len(Q().subset_by_overlaps(S())) == 3
    assert len(Q().subset_by_overlaps(S(), invert=True)) == 2


def test_overlap_types_and_maxgap_use_the_iranges_predicate():
    q = GRanges("chr1", IRanges(1, 10), "+")
    s = GRanges("chr1", IRanges([5, 11, 20], width=3), "+")
    assert pairs(q.find_overlaps(s, type="within")) == []
    assert pairs(s.find_overlaps(q, type="within")) == [(0, 0)]
    assert pairs(q.find_overlaps(s, maxgap=0)) == [(0, 0), (0, 1)]  # adjacent 11-13
    assert pairs(q.find_overlaps(s, maxgap=10)) == [(0, 0), (0, 1), (0, 2)]
    assert pairs(q.find_overlaps(s, minoverlap=4)) == []
    assert pairs(q.find_overlaps(s, minoverlap=3)) == [(0, 0)]
    with pytest.raises(ValueError, match="default value"):
        q.find_overlaps(s, maxgap=0, minoverlap=2)


def test_self_overlaps():
    # R: findOverlaps(q, drop.self=TRUE) -> (1,3) (2,3) (3,1) (3,2)
    h = Q().find_overlaps(drop_self=True)
    assert pairs(h) == [(0, 2), (1, 2), (2, 0), (2, 1)]
    assert Q().find_overlaps(select="first").tolist() == [0, 1, 0, 3, 4]


def test_different_seqlevels_are_matched_by_name():
    q = GRanges("chrB", IRanges(1, 10), seqinfo=Seqinfo(["chrA", "chrB"]))
    s = GRanges(["chrB", "chrA"], IRanges(5, 6), seqinfo=Seqinfo(["chrB", "chrA"]))
    assert pairs(q.find_overlaps(s)) == [(0, 0)]
    other = GRanges("chrX", IRanges(1, 10))
    assert q.count_overlaps(other).tolist() == [0]


def test_empty_inputs():
    e = GRanges.empty()
    assert len(Q().find_overlaps(e)) == 0 and len(e.find_overlaps(S())) == 0
    assert Q().find_overlaps(e, select="first").tolist() == [NA] * 5


# ---- nearest family -----------------------------------------------------------------
def NQ():
    return GRanges(["chr1", "chr1", "chr2", "chr1"], IRanges([10, 10, 10, 100], width=5), ["+", "-", "+", "*"])


def NS():
    return GRanges(["chr1", "chr1", "chr1", "chr2", "chr1"], IRanges([1, 30, 50, 40, 200], width=5),
                   ["+", "+", "-", "+", "*"])


def test_precede_follow_nearest_are_strand_aware():
    # R (1-based): precede 2 NA 4 3 ; follow 1 3 NA 2 ; nearest 1 3 4 3
    assert NQ().precede(NS()).tolist() == [1, NA, 3, 2]
    assert NQ().follow(NS()).tolist() == [0, 2, NA, 1]
    assert NQ().nearest(NS()).tolist() == [0, 2, 3, 2]


def test_precede_follow_nearest_ignore_strand():
    # R: precede 2 2 4 5 ; follow 1 1 NA 3 ; nearest 1 1 4 3
    assert NQ().precede(NS(), ignore_strand=True).tolist() == [1, 1, 3, 4]
    assert NQ().follow(NS(), ignore_strand=True).tolist() == [0, 0, NA, 2]
    assert NQ().nearest(NS(), ignore_strand=True).tolist() == [0, 0, 3, 2]


def test_star_query_takes_the_closer_direction():
    # '*' query at 100-104: + subjects to its right (200), - subjects to its left (50-54, closer)
    x = GRanges("chr1", IRanges(100, 104), "*")
    s = GRanges("chr1", IRanges([50, 200], width=5), ["-", "+"])
    assert x.precede(s).tolist() == [0]  # R: 1 ... the - subject on the left wins
    # R: gaps of 45 (the - subject at 50-54) and 45 (the + subject at 150): tie keeps both, "first" -> lowest index
    s2 = GRanges("chr1", IRanges([50, 150], width=5), ["-", "+"])
    assert x.precede(s2).tolist() == [0]
    assert pairs(x.precede(s2, select="all")) == [(0, 0), (0, 1)]
    # the + subject at 145-150 is closer (gap 40) and wins
    s3 = GRanges("chr1", IRanges([50, 145], width=[5, 6]), ["-", "+"])
    assert x.precede(s3).tolist() == [1]
    assert pairs(x.precede(s3, select="all")) == [(0, 1)]


def test_precede_ties_first_follow_ties_last():
    s = GRanges("chr1", IRanges([30, 30, 5, 5], width=5), "+")
    x = GRanges("chr1", IRanges(10, 14), "+")
    assert x.precede(s).tolist() == [0]
    assert pairs(x.precede(s, select="all")) == [(0, 0), (0, 1)]
    y = GRanges("chr1", IRanges(100, 104), "+")
    assert y.follow(s).tolist() == [1]
    assert pairs(y.follow(s, select="all")) == [(0, 0), (0, 1)]


def test_precede_follow_self_and_all():
    x = GRanges("chr1", IRanges([1, 5, 10], [3, 6, 12]), "+")
    assert x.precede().tolist() == [1, 2, NA]
    assert x.follow().tolist() == [NA, 0, 1]
    assert x.nearest().tolist() == [1, 0, 1]  # R: 2 1 2
    assert pairs(x.nearest(select="all")) == [(0, 1), (1, 0), (2, 1)]
    assert pairs(x.precede(select="all")) == [(0, 1), (1, 2)]


def test_nearest_overlap_or_adjacent_wins():
    x = GRanges("chr1", IRanges(10, 12), "+")
    s = GRanges("chr1", IRanges(13, 15), "-")
    # R: strand-incompatible neighbours don't count -> NA ; ignoring strand it's adjacent
    assert x.nearest(s).tolist() == [NA]
    assert x.nearest(s, ignore_strand=True).tolist() == [0]


def test_nearest_equidistant():
    x = GRanges("chr1", IRanges(5, 6), "+")
    # R: arbitrary -> lower subject index among the two sides; all -> precede side only
    s = GRanges("chr1", IRanges([8, 1], [9, 3]), "+")
    assert x.nearest(s).tolist() == [0]
    assert pairs(x.nearest(s, select="all")) == [(0, 0)]
    s = GRanges("chr1", IRanges([1, 8], [3, 9]), "+")
    assert x.nearest(s).tolist() == [0]
    assert pairs(x.nearest(s, select="all")) == [(0, 1)]


def test_nearest_empty_subject():
    x = GRanges("chr1", IRanges(5, 6), "+")
    assert x.nearest(GRanges.empty()).tolist() == [NA]
    assert x.precede(GRanges.empty()).tolist() == [NA]


def test_distance_na_rules():
    # R: distance(q[1:4], s[1:4]) -> 4 NA NA NA
    assert NQ().distance(NS()[:4]).tolist() == [4, NA, NA, NA]
    # R: ignore.strand, recycled length-1 subject -> 4 4 NA 94
    assert NQ().distance(NS()[:1], ignore_strand=True).tolist() == [4, 4, NA, 94]
    # '*' is compatible with both strands
    a = GRanges("chr1", IRanges(1, 5), "*")
    b = GRanges("chr1", IRanges(10, 12), ["+", "-"][:1])
    assert a.distance(b).tolist() == [4]
    assert len(GRanges.empty().distance(b)) == 0


# ---- coverage ----------------------------------------------------------------------------
def test_coverage_per_seqlevel():
    si = Seqinfo(["chr1", "chr2", "chr3"], [20, 10, None])
    g = GRanges(["chr1", "chr1", "chr2", "chr3"], IRanges([2, 5, 3, 4], width=[5, 3, 2, 6]),
                ["+", "-", "+", "*"], seqinfo=si)
    cv = g.coverage()
    # R: chr1 lengths 1 3 2 1 13 / values 0 1 2 1 0 ; chr2 20? no: seqlength 10 ; chr3 unknown -> max end 9
    assert list(cv) == ["chr1", "chr2", "chr3"]
    assert cv["chr1"].lengths.tolist() == [1, 3, 2, 1, 13] and cv["chr1"].values.tolist() == [0, 1, 2, 1, 0]
    assert len(cv["chr2"]) == 10 and cv["chr2"].lengths.tolist() == [2, 2, 6]
    assert len(cv["chr3"]) == 9 and cv["chr3"].lengths.tolist() == [3, 6]


def test_coverage_width_shift_unused_levels():
    si = Seqinfo(["chr1", "chr2", "chr3"], [20, 10, None])
    g = GRanges(["chr1", "chr1", "chr2", "chr3"], IRanges([2, 5, 3, 4], width=[5, 3, 2, 6]), seqinfo=si)
    assert {k: len(v) for k, v in g.coverage(width=30).items()} == {"chr1": 30, "chr2": 30, "chr3": 30}
    cs = g.coverage(shift=1)
    assert cs["chr1"].lengths.tolist() == [2, 3, 2, 1, 12]  # R
    assert len(cs["chr3"]) == 10
    h = GRanges(["chr1", "chr2"], IRanges([2, 3], width=[5, 2]), seqinfo=Seqinfo(["chr1", "chr2", "chr9"]))
    cv = h.coverage()
    assert len(cv["chr9"]) == 0 and len(cv["chr1"]) == 6 and len(cv["chr2"]) == 4  # R


def test_empty_range_is_its_own_neighbour_like_r():
    # R: g <- GRanges("chr1", IRanges(5,4), "+"); nearest(g), precede(g), follow(g) -> 1 (review finding 1)
    g = GRanges("chr1", IRanges(5, 4), "+")
    assert g.nearest().tolist() == [0]
    assert g.precede().tolist() == [0] and g.follow().tolist() == [0]


def test_distance_rejects_incompatible_lengths():
    # review finding 2: only length-1 operands may be recycled
    a = GRanges("chr1", IRanges([1, 100], width=1))
    b = GRanges("chr1", IRanges([10, 20, 30], width=1))
    with pytest.raises(ValueError, match="same length"):
        a.distance(b)
    with pytest.raises(ValueError, match="same length"):
        b.distance(a)
    assert a[:1].distance(b).tolist() == [8, 18, 28]
    assert b.distance(a[:1]).tolist() == [8, 18, 28]
