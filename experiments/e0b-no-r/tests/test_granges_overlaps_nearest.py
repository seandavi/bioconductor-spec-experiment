
import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from ranges_e0b import NO_HIT, GRanges, IRanges, Seqinfo


def g(seq: list[str], s: list[int], e: list[int], strand: list[str]) -> GRanges:
    return GRanges(seq, IRanges(start=s, end=e), strand)


def test_overlaps_same_seqlevel_only() -> None:
    q = g(["a", "b"], [1, 1], [10, 10], ["*", "*"])
    s = g(["b", "a"], [5, 5], [6, 6], ["*", "*"])
    assert q.find_overlaps(s).pairs() == [(0, 1), (1, 0)]


def test_strand_compatibility_matrix() -> None:
    q = g(["a"] * 3, [1, 1, 1], [10, 10, 10], ["+", "-", "*"])
    s = g(["a"] * 3, [5, 5, 5], [6, 6, 6], ["+", "-", "*"])
    # + overlaps + and *, - overlaps - and *, * overlaps all
    assert q.find_overlaps(s).pairs() == [
        (0, 0), (0, 2), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)
    ]
    assert len(q.find_overlaps(s, ignore_strand=True)) == 9


def test_select_count_any_subset() -> None:
    q = g(["a"] * 3, [1, 1, 50], [10, 10, 60], ["+", "-", "+"])
    s = g(["a"] * 2, [5, 5], [6, 8], ["-", "-"])
    assert q.find_overlaps(s, select="first").tolist() == [NO_HIT, 0, NO_HIT]
    assert q.find_overlaps(s, select="last").tolist() == [NO_HIT, 1, NO_HIT]
    assert q.count_overlaps(s).tolist() == [0, 2, 0]
    assert q.count_overlaps(s, ignore_strand=True).tolist() == [2, 2, 0]
    assert q.overlaps_any(s).tolist() == [False, True, False]
    assert q.subset_by_overlaps(s).strand == ["-"]
    assert q.subset_by_overlaps(s, invert=True).start.tolist() == [1, 50]


def test_start_end_flip_on_minus() -> None:
    q = g(["a", "a"], [1, 1], [10, 10], ["+", "-"])
    s = g(["a"], [1], [20], ["*"])
    # start type: + compares starts (equal) ; - compares ends (10 vs 20: no)
    assert q.find_overlaps(s, type="start").pairs() == [(0, 0)]
    # end type: + compares ends (no); - compares starts (equal)
    assert q.find_overlaps(s, type="end").pairs() == [(1, 0)]
    assert q.find_overlaps(s, type="end", ignore_strand=True).pairs() == []
    assert q.find_overlaps(s, type="within").pairs() == [(0, 0), (1, 0)]
    assert q.find_overlaps(s, type="equal").pairs() == []


def test_overlaps_merges_seqinfo_by_name() -> None:
    q = GRanges(["a"], IRanges(start=1, end=10), seqinfo=Seqinfo(["a", "b"]))
    s = GRanges(["a"], IRanges(start=5, end=6), seqinfo=Seqinfo(["b", "a"]))
    assert q.find_overlaps(s).pairs() == [(0, 0)]


def test_overlaps_maxgap_minoverlap() -> None:
    q = g(["a"], [1], [5], ["+"])
    s = g(["a", "a"], [6, 4], [9, 9], ["+", "+"])
    assert q.find_overlaps(s, maxgap=0).pairs() == [(0, 0), (0, 1)]
    assert q.find_overlaps(s, minoverlap=2).pairs() == [(0, 1)]


# --- nearest ------------------------------------------------------------------


def test_precede_follow_strand_direction() -> None:
    s = g(["a", "a"], [10, 30], [12, 32], ["*", "*"])
    q = g(["a", "a", "a"], [20, 20, 20], [22, 22, 22], ["+", "-", "*"])
    # + looks to higher coordinates; - to lower; * treated as +
    assert q.precede(s).tolist() == [1, 0, 1]
    assert q.follow(s).tolist() == [0, 1, 0]


def test_precede_requires_same_seqlevel_and_compatible_strand() -> None:
    s = g(["a", "b", "a"], [30, 30, 40], [31, 31, 41], ["+", "+", "-"])
    q = g(["a"], [1], [2], ["+"])
    assert q.precede(s).tolist() == [0]
    q2 = g(["a"], [1], [2], ["-"])
    assert q2.precede(s).tolist() == [NO_HIT]  # lower coordinates only: nothing
    assert q2.follow(s).tolist() == [2]  # '-' looks higher; s[0] is '+' (incompatible), s[2] is '-'
    assert q2.follow(s, ignore_strand=True).tolist() == [NO_HIT]  # strand ignored -> treated as '+'
    assert q2.precede(s, ignore_strand=True).tolist() == [0]
    assert q.precede(s, ignore_strand=True).tolist() == [0]


def test_follow_minus_strand_picks_higher() -> None:
    s = g(["a", "a"], [30, 50], [31, 51], ["-", "-"])
    q = g(["a"], [1], [2], ["-"])
    assert q.follow(s).tolist() == [0]
    assert q.precede(s).tolist() == [NO_HIT]


def test_nearest_per_seqlevel_and_strand() -> None:
    s = g(["a", "a", "b"], [100, 10, 3], [110, 12, 4], ["+", "-", "+"])
    q = g(["a", "a", "b", "c"], [20, 20, 1, 1], [22, 22, 2, 2], ["+", "-", "*", "+"])
    # + query: only '+'/'*' subjects -> idx 0 ; - query -> idx 1 ; b -> idx 2 ; c -> none
    assert q.nearest(s).tolist() == [0, 1, 2, NO_HIT]
    assert q.nearest(s, ignore_strand=True).tolist() == [1, 1, 2, NO_HIT]
    assert q.nearest(s, select="all").pairs() == [(0, 0), (1, 1), (2, 2)]


def test_nearest_self() -> None:
    x = g(["a", "a", "a"], [1, 10, 30], [5, 12, 40], ["+", "+", "-"])
    assert x.nearest().tolist() == [1, 0, NO_HIT]
    assert x.nearest(ignore_strand=True).tolist() == [1, 0, 1]
    assert x.precede().tolist() == [1, NO_HIT, NO_HIT]
    assert x.follow().tolist() == [NO_HIT, 0, NO_HIT]


def test_distance() -> None:
    a = g(["a", "a", "a", "b"], [1, 1, 1, 1], [5, 5, 5, 5], ["+", "+", "+", "*"])
    b = g(["a", "a", "b", "b"], [10, 3, 10, 10], [12, 4, 12, 12], ["+", "+", "+", "-"])
    assert a.distance(b).tolist() == [5, 0, NO_HIT, 5]
    c = g(["a"], [10], [12], ["-"])
    assert a[[0]].distance(c).tolist() == [NO_HIT]
    assert a[[0]].distance(c, ignore_strand=True).tolist() == [5]
    with pytest.raises(ValueError):
        a.distance(c)


# --- coverage -----------------------------------------------------------------


def test_coverage_dict_over_seqlevels() -> None:
    x = GRanges(["a", "a", "b"], IRanges(start=[1, 3, 2], end=[4, 5, 2]), ["+", "-", "+"],
                seqinfo=Seqinfo(["a", "b", "c"], [6, 4, 3]))
    cov = x.coverage()
    assert list(cov) == ["a", "b", "c"]
    assert cov["a"].to_numpy().tolist() == [1, 1, 2, 2, 1, 0]  # strand ignored, seqlength 6
    assert cov["b"].to_numpy().tolist() == [0, 1, 0, 0]
    assert cov["c"].to_numpy().tolist() == [0, 0, 0]
    assert x.coverage(width=5)["b"].to_numpy().tolist() == [0, 1, 0, 0, 0]
    assert x.coverage(width={"a": 2})["a"].to_numpy().tolist() == [1, 1]
    assert x.coverage(weight=[1, 2, 3])["a"].to_numpy().tolist() == [1, 1, 3, 3, 2, 0]
    with pytest.raises(ValueError):
        x.coverage(weight=[1])


def test_coverage_unknown_length_uses_max_end() -> None:
    x = GRanges(["a"], IRanges(start=[2], end=[4]))
    assert x.coverage()["a"].to_numpy().tolist() == [0, 1, 1, 1]


# --- property: GRanges == IRanges per (seq, strand) --------------------------------

rec = st.tuples(st.sampled_from(["a", "b"]), st.sampled_from("+-*"), st.integers(1, 25), st.integers(1, 6))


def build(rs: list[tuple[str, str, int, int]]) -> GRanges:
    if not rs:
        return GRanges([], IRanges(), seqinfo=Seqinfo(["a", "b"]))
    return GRanges([r[0] for r in rs], IRanges(start=[r[2] for r in rs], width=[r[3] for r in rs]),
                   [r[1] for r in rs], seqinfo=Seqinfo(["a", "b"]))


@given(st.lists(rec, max_size=12), st.lists(rec, max_size=12), st.booleans(),
       st.sampled_from(["any", "within", "equal"]), st.integers(-1, 3))
def test_overlaps_match_pairwise_brute_force(xs: list, ys: list, ignore: bool, type: str, maxgap: int) -> None:  # type: ignore[type-arg]
    got = build(xs).find_overlaps(build(ys), type=type, maxgap=maxgap, ignore_strand=ignore).pairs()
    expect = []
    for i, (qa, qst, qs, qw) in enumerate(xs):
        for j, (sa, sst, ss, sw) in enumerate(ys):
            if qa != sa or (not ignore and not (qst == sst or "*" in (qst, sst))):
                continue
            m = max(maxgap, 0)
            qe, se = qs + qw - 1, ss + sw - 1
            ok = {
                "any": max(qs, ss) - min(qe, se) - 1 <= maxgap,
                "within": qs >= ss and qe <= se and (maxgap < 0 or (qs - ss) + (se - qe) <= maxgap),
                "equal": abs(qs - ss) <= m and abs(qe - se) <= m,
            }[type]
            if ok:
                expect.append((i, j))
    assert got == expect


@given(st.lists(rec, min_size=1, max_size=10), st.lists(rec, max_size=10), st.booleans())
def test_nearest_matches_brute_force(xs: list, ys: list, ignore: bool) -> None:  # type: ignore[type-arg]
    got = build(xs).nearest(build(ys), select="all", ignore_strand=ignore).pairs()
    expect = []
    for i, (qa, qst, qs, qw) in enumerate(xs):
        best, ids = None, []
        for j, (sa, sst, ss, sw) in enumerate(ys):
            if qa != sa or (not ignore and not (qst == sst or "*" in (qst, sst))):
                continue
            d = max(max(qs, ss) - min(qs + qw - 1, ss + sw - 1) - 1, 0)
            if best is None or d < best:
                best, ids = d, [j]
            elif d == best:
                ids.append(j)
        expect += [(i, j) for j in ids]
    assert got == expect


@given(st.lists(rec, max_size=12), st.lists(rec, max_size=12))
def test_precede_follow_strand_brute_force(xs: list, ys: list) -> None:  # type: ignore[type-arg]
    q, s = build(xs), build(ys)
    pre, fol = q.precede(s), q.follow(s)
    for i, (qa, qst, qs, qw) in enumerate(xs):
        qe = qs + qw - 1
        cands_pre, cands_fol = [], []
        for j, (sa, sst, ss, sw) in enumerate(ys):
            if qa != sa or not (qst == sst or "*" in (qst, sst)):
                continue
            se = ss + sw - 1
            minus = qst == "-"
            if (ss > qe and not minus) or (se < qs and minus):
                cands_pre.append((ss if not minus else -se, j))
            if (se < qs and not minus) or (ss > qe and minus):
                cands_fol.append((-se if not minus else ss, j))
        assert pre[i] == (min(cands_pre)[1] if cands_pre else NO_HIT)
        assert fol[i] == (min(cands_fol)[1] if cands_fol else NO_HIT)


@given(st.lists(rec, max_size=12))
def test_coverage_matches_dense(xs: list) -> None:  # type: ignore[type-arg]
    cov = build(xs).coverage()
    for name in ("a", "b"):
        dense = np.zeros(max([r[2] + r[3] - 1 for r in xs if r[0] == name] or [0]), dtype=np.int64)
        for a, _, s, w in xs:
            if a == name:
                dense[s - 1 : s - 1 + w] += 1
        assert cov[name].to_numpy().tolist() == dense.tolist()
