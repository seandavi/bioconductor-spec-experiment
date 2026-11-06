import itertools

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from ranges_e0 import NO_MATCH, Hits, IRanges, _overlaps


def hit_starts(subject, width, **kw):
    """Query starts 1..12 of the given width that hit `subject` (hand-checked in R)."""
    q = IRanges(np.arange(1, 13), width=width)
    return [int(s) for s, c in zip(q.start, q.count_overlaps(subject, **kw)) if c > 0]


# ---- type="any" -------------------------------------------------------------
def test_any_default_closed_interval_overlap():
    s = IRanges(4, 6)
    # R: subject 4-6; width-0 queries starting at 5,6 overlap (inside), at 4 and 7 do not
    assert hit_starts(s, 0) == [5, 6]
    assert hit_starts(s, 1) == [4, 5, 6]
    assert hit_starts(s, 3) == [2, 3, 4, 5, 6]


def test_any_empty_subject():
    s = IRanges(5, 4)
    assert hit_starts(s, 0) == []
    assert hit_starts(s, 2) == [4]  # query 4-5 contains the point between 4 and 5
    assert hit_starts(s, 3) == [3, 4]


def test_any_maxgap():
    s = IRanges(4, 6)
    # R: maxgap=0 includes adjacent ranges, maxgap=1 one more position
    assert hit_starts(s, 0, maxgap=0) == [4, 5, 6, 7]
    assert hit_starts(s, 1, maxgap=0) == [3, 4, 5, 6, 7]
    assert hit_starts(s, 0, maxgap=1) == [3, 4, 5, 6, 7, 8]
    assert hit_starts(s, 3, maxgap=1) == [1, 2, 3, 4, 5, 6, 7, 8]
    assert hit_starts(IRanges(5, 4), 0, maxgap=1) == [4, 5, 6]


def test_any_minoverlap():
    s = IRanges(4, 6)
    assert hit_starts(s, 0, minoverlap=1) == []
    assert hit_starts(s, 1, minoverlap=1) == [4, 5, 6]
    assert hit_starts(s, 2, minoverlap=2) == [4, 5]
    assert hit_starts(s, 3, minoverlap=2) == [3, 4, 5]


def test_any_maxgap_and_minoverlap_together_is_an_error():
    # R: 'at least one of maxgap and minoverlap must be set to its default value'
    with pytest.raises(ValueError, match="default value"):
        IRanges(1, 3).find_overlaps(IRanges(2, 4), maxgap=0, minoverlap=2)
    with pytest.raises(ValueError):
        IRanges(1, 3).find_overlaps(IRanges(2, 4), maxgap=1, minoverlap=2)
    IRanges(1, 3).find_overlaps(IRanges(2, 4), maxgap=-1, minoverlap=2)  # fine


def test_argument_validation():
    q = IRanges(1, 3)
    with pytest.raises(ValueError, match="maxgap"):
        q.find_overlaps(q, maxgap=-2)
    with pytest.raises(ValueError, match="minoverlap"):
        q.find_overlaps(q, minoverlap=-1)
    with pytest.raises(ValueError, match="type"):
        q.find_overlaps(q, type="foo")
    with pytest.raises(ValueError, match="select"):
        q.find_overlaps(q, select="foo")


# ---- special types (subject 4-9) -------------------------------------------
S = IRanges(4, 9)


def test_type_start():
    # R: start match: query start must equal 4, any width (incl. width 0)
    assert hit_starts(S, 0, type="start") == [4]
    assert hit_starts(S, 7, type="start") == [4]
    # maxgap=2: starts within 2 of 4, any width
    assert hit_starts(S, 3, type="start", maxgap=2) == [2, 3, 4, 5, 6]
    assert hit_starts(S, 0, type="start", maxgap=2) == [2, 3, 4, 5, 6]


def test_type_end():
    # width 0 at start 10 is (10,9): ends at 9
    assert hit_starts(S, 0, type="end") == [10]
    assert hit_starts(S, 5, type="end") == [5]
    assert hit_starts(S, 7, type="end") == [3]
    assert hit_starts(S, 0, type="end", maxgap=2) == [8, 9, 10, 11, 12]


def test_type_equal():
    assert hit_starts(S, 6, type="equal") == [4]
    assert hit_starts(S, 5, type="equal") == []
    # maxgap=1: both ends within 1
    assert hit_starts(S, 4, type="equal", maxgap=1) == [5]
    assert hit_starts(S, 6, type="equal", maxgap=1) == [3, 4, 5]
    assert hit_starts(S, 8, type="equal", maxgap=1) == [3]
    assert hit_starts(S, 4, type="equal", maxgap=1, minoverlap=5) == []


def test_type_within():
    # R: query wholly inside subject 4-9; zero-width at 4..10 counts as within
    assert hit_starts(S, 0, type="within") == [4, 5, 6, 7, 8, 9, 10]
    assert hit_starts(S, 3, type="within") == [4, 5, 6, 7]
    assert hit_starts(S, 6, type="within") == [4]
    assert hit_starts(S, 7, type="within") == []


def test_type_within_maxgap_limits_subject_excess_but_zero_means_unlimited():
    # R quirk: maxgap 0 behaves like -1 (no limit); maxgap>0 bounds sw - qw
    assert hit_starts(S, 0, type="within", maxgap=0) == [4, 5, 6, 7, 8, 9, 10]
    assert hit_starts(S, 3, type="within", maxgap=0) == [4, 5, 6, 7]
    assert hit_starts(S, 4, type="within", maxgap=1) == []
    assert hit_starts(S, 5, type="within", maxgap=1) == [4, 5]
    assert hit_starts(S, 4, type="within", maxgap=2) == [4, 5, 6]
    assert hit_starts(S, 3, type="within", maxgap=2) == []


def test_type_within_minoverlap():
    assert hit_starts(S, 2, type="within", minoverlap=3) == []
    assert hit_starts(S, 3, type="within", minoverlap=3) == [4, 5, 6, 7]


def test_start_with_minoverlap_and_maxgap():
    assert hit_starts(S, 2, type="start", minoverlap=3) == []
    assert hit_starts(S, 3, type="start", minoverlap=3) == [4]
    # R: starts within 2 of 4 and overlap length >= 3
    assert hit_starts(S, 3, type="start", maxgap=2, minoverlap=3) == [4, 5, 6]
    assert hit_starts(S, 5, type="start", maxgap=2, minoverlap=3) == [2, 3, 4, 5, 6]


def test_empty_subject_special_types():
    e = IRanges(6, 5)
    assert hit_starts(e, 0, type="start") == [6] and hit_starts(e, 3, type="start") == [6]
    assert hit_starts(e, 0, type="within") == [6]
    assert hit_starts(e, 1, type="within") == []


# ---- select / Hits / helpers -----------------------------------------------
Q = IRanges([5, 1, 3], [9, 4, 3])
SUB = IRanges([8, 2, 1, 3], [12, 3, 2, 6])


def test_find_overlaps_all_hits():
    # R (1-based): (1,1) (1,4) (2,2) (2,3) (2,4) (3,2) (3,4)
    h = Q.find_overlaps(SUB)
    assert isinstance(h, Hits)
    assert h.to_numpy().tolist() == [[0, 0], [0, 3], [1, 1], [1, 2], [1, 3], [2, 1], [2, 3]]
    assert (h.query_length, h.subject_length) == (3, 4)


def test_find_overlaps_select():
    assert Q.find_overlaps(SUB, select="first").tolist() == [0, 1, 1]
    assert Q.find_overlaps(SUB, select="last").tolist() == [3, 3, 3]
    arb = Q.find_overlaps(SUB, select="arbitrary")
    hits = {tuple(p) for p in Q.find_overlaps(SUB).to_numpy().tolist()}
    assert all((i, int(j)) in hits for i, j in enumerate(arb))
    assert Q.find_overlaps(IRanges.empty(), select="first").tolist() == [NO_MATCH] * 3
    assert len(Q.find_overlaps(IRanges.empty())) == 0


def test_within_hits():
    # R: findOverlaps(q, s, type="within") -> (3,2) (3,4)
    assert Q.find_overlaps(SUB, type="within").to_numpy().tolist() == [[2, 1], [2, 3]]
    assert Q.find_overlaps(SUB, type="start", maxgap=0).to_numpy().tolist() == [[1, 2], [2, 3]]


def test_count_any_subset():
    assert Q.count_overlaps(SUB).tolist() == [2, 3, 2]
    assert Q.overlaps_any(SUB).tolist() == [True, True, True]
    assert Q.overlaps_any(IRanges(100, 101)).tolist() == [False] * 3
    assert len(Q.subset_by_overlaps(SUB, invert=True)) == 0
    assert Q.subset_by_overlaps(IRanges(8, 9)).to_tuples() == [(5, 9)]
    assert Q.subset_by_overlaps(IRanges(8, 9), invert=True).to_tuples() == [(1, 4), (3, 3)]


def test_subset_by_overlaps_keeps_names_and_mcols():
    q = Q.with_names(["a", "b", "c"]).with_mcols({"v": [1, 2, 3]})
    r = q.subset_by_overlaps(IRanges(8, 9))
    assert r.names == ["a"] and r.mcols["v"].to_list() == [1]


def test_self_overlaps():
    x = IRanges([1, 2, 2, 6, 3], [5, 3, 4, 7, 3])
    pairs = lambda h: sorted(map(tuple, h.to_numpy().tolist()))  # noqa: E731
    full = pairs(x.find_overlaps())
    assert len(full) == 17  # R
    assert (0, 0) in full and (3, 3) in full
    assert len(pairs(x.find_overlaps(drop_self=True))) == 12
    # R: drop.redundant keeps one of A->B / B->A (and the self hits): 11 pairs
    red = pairs(x.find_overlaps(drop_redundant=True))
    assert len(red) == 11 and all(i <= j for i, j in red)
    both = pairs(x.find_overlaps(drop_self=True, drop_redundant=True))
    assert len(both) == 6 and all(i < j for i, j in both)
    assert x.find_overlaps(select="first").tolist() == [0, 0, 0, 3, 0]
    with pytest.raises(ValueError):
        x.find_overlaps(x, drop_self=True)


def test_hits_are_sorted_by_query_then_subject():
    h = IRanges([1, 1, 1], [5, 5, 5]).find_overlaps(IRanges([3, 1, 2], [4, 6, 3]))
    assert h.to_numpy().tolist() == [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]


# ---- kernel vs brute force --------------------------------------------------
def brute(q: IRanges, s: IRanges, type: str, maxgap: int, minoverlap: int) -> list[tuple[int, int]]:
    """Independent O(n*m) reading of the documented predicate, with explicit loops."""
    out = []
    gg = max(maxgap, 0)
    for i, (qs, qe) in enumerate(q.to_tuples()):
        for j, (ss, se) in enumerate(s.to_tuples()):
            if type == "any":
                ok = qs <= se + maxgap + 1 and ss <= qe + maxgap + 1
            elif type == "start":
                ok = abs(qs - ss) <= gg
            elif type == "end":
                ok = abs(qe - se) <= gg
            elif type == "equal":
                ok = abs(qs - ss) <= gg and abs(qe - se) <= gg
            else:
                ok = ss <= qs and qe <= se
                if ok and maxgap > 0:
                    ok = (se - ss + 1) - (qe - qs + 1) <= maxgap
            if ok and minoverlap > 0:
                ok = min(qe, se) - max(qs, ss) + 1 >= minoverlap
            if ok:
                out.append((i, j))
    return out


ranges_st = st.lists(st.tuples(st.integers(-5, 25), st.integers(0, 8)), min_size=0, max_size=12).map(
    lambda xs: IRanges([a for a, _ in xs], width=[w for _, w in xs]) if xs else IRanges.empty()
)


@settings(max_examples=300, deadline=None)
@given(
    q=ranges_st,
    s=ranges_st,
    type=st.sampled_from(_overlaps.TYPES),
    maxgap=st.integers(-1, 4),
    minoverlap=st.integers(0, 4),
)
def test_kernel_matches_bruteforce(q, s, type, maxgap, minoverlap):
    if type == "any" and maxgap != -1 and minoverlap != 0:
        return
    h = q.find_overlaps(s, type=type, maxgap=maxgap, minoverlap=minoverlap)
    assert list(h) == sorted(brute(q, s, type, maxgap, minoverlap))
    counts = q.count_overlaps(s, type=type, maxgap=maxgap, minoverlap=minoverlap)
    assert counts.tolist() == np.bincount(h.query_hits, minlength=len(q)).tolist()


def test_chunked_path_matches_unchunked(monkeypatch):
    rng = np.random.default_rng(0)
    q = IRanges(rng.integers(1, 200, 60), width=rng.integers(0, 15, 60))
    s = IRanges(rng.integers(1, 200, 50), width=rng.integers(0, 15, 50))
    expected = list(q.find_overlaps(s))
    monkeypatch.setattr(_overlaps, "_CHUNK", 7)
    assert list(q.find_overlaps(s)) == expected
