import numpy as np
import pytest
from hypothesis import given, strategies as st

from helpers import iranges_st
from ranges_e0b import NO_HIT, Hits, IRanges

S = IRanges(start=[1, 6, 11], end=[5, 10, 15])
Q = IRanges(start=[4, 5, 20], end=[6, 5, 25])


def test_any() -> None:
    # R: findOverlaps(Q, S)  (1-based: 1-1, 1-2, 2-1)
    h = Q.find_overlaps(S)
    assert isinstance(h, Hits)
    assert h.pairs() == [(0, 0), (0, 1), (1, 0)]
    assert (h.query_length, h.subject_length) == (3, 3)


def test_maxgap() -> None:
    # R: findOverlaps(Q, S, maxgap=0)  -- adjacency counts
    assert Q.find_overlaps(S, maxgap=0).pairs() == [(0, 0), (0, 1), (1, 0), (1, 1)]
    assert Q.find_overlaps(S, maxgap=4).pairs() == [
        (0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 2)
    ]


def test_minoverlap() -> None:
    # R: findOverlaps(Q, S, minoverlap=2)
    assert Q.find_overlaps(S, minoverlap=2).pairs() == [(0, 0)]
    assert Q.find_overlaps(S, minoverlap=1).pairs() == Q.find_overlaps(S).pairs()


def test_any_maxgap_and_minoverlap_conflict() -> None:
    with pytest.raises(ValueError):
        Q.find_overlaps(S, maxgap=1, minoverlap=1)
    with pytest.raises(ValueError):
        Q.find_overlaps(S, maxgap=-2)
    with pytest.raises(ValueError):
        Q.find_overlaps(S, type="bogus")


def test_type_start_end_within_equal() -> None:
    assert Q.find_overlaps(S, type="start").pairs() == []
    assert Q.find_overlaps(S, type="start", maxgap=1).pairs() == [(1, 1)]
    assert Q.find_overlaps(S, type="end").pairs() == [(1, 0)]
    assert Q.find_overlaps(S, type="within").pairs() == [(1, 0)]
    q = IRanges(start=[6, 6], end=[10, 9])
    assert q.find_overlaps(S, type="equal").pairs() == [(0, 1)]
    assert q.find_overlaps(S, type="equal", maxgap=1).pairs() == [(0, 1), (1, 1)]
    assert q.find_overlaps(S, type="start").pairs() == [(0, 1), (1, 1)]


def test_within_nested_and_minoverlap() -> None:
    s = IRanges(start=[1, 3], end=[10, 4])
    q = IRanges(start=[3, 2], end=[4, 4])
    assert q.find_overlaps(s, type="within").pairs() == [(0, 0), (0, 1), (1, 0)]
    assert q.find_overlaps(s, type="within", minoverlap=2).pairs() == [(0, 0), (0, 1), (1, 0)]
    assert q.find_overlaps(s, type="within", minoverlap=3).pairs() == [(1, 0)]


def test_select() -> None:
    assert Q.find_overlaps(S, select="first").tolist() == [0, 0, NO_HIT]
    assert Q.find_overlaps(S, select="last").tolist() == [1, 0, NO_HIT]
    assert Q.find_overlaps(S, select="arbitrary").tolist() == [0, 0, NO_HIT]
    with pytest.raises(ValueError):
        Q.find_overlaps(S, select="nope")


def test_count_any_subset() -> None:
    assert Q.count_overlaps(S).tolist() == [2, 1, 0]
    assert Q.overlaps_any(S).tolist() == [True, True, False]
    assert Q.subset_by_overlaps(S).start.tolist() == [4, 5]
    assert Q.subset_by_overlaps(S, invert=True).start.tolist() == [20]


def test_subset_by_overlaps_keeps_names() -> None:
    q = Q.with_names(["a", "b", "c"])
    assert q.subset_by_overlaps(S).names == ["a", "b"]


def test_self_overlaps_default_subject() -> None:
    x = IRanges(start=[1, 3, 10], end=[5, 8, 12])
    assert x.find_overlaps().pairs() == [(0, 0), (0, 1), (1, 0), (1, 1), (2, 2)]


def test_empty_inputs() -> None:
    assert len(IRanges().find_overlaps(S)) == 0
    assert len(Q.find_overlaps(IRanges())) == 0
    assert Q.find_overlaps(IRanges(), select="first").tolist() == [NO_HIT] * 3
    assert Q.count_overlaps(IRanges()).tolist() == [0, 0, 0]


def test_zero_width_query() -> None:
    # UNVERIFIED: empty range at p overlaps [a,b] iff a < p <= b
    s = IRanges(start=1, end=10)
    q = IRanges(start=[5, 1, 10, 11], end=[4, 0, 9, 10])
    assert q.overlaps_any(s).tolist() == [True, False, True, False]
    assert q.overlaps_any(s, minoverlap=1).tolist() == [False] * 4
    # touching ends hit once gaps up to 0 are allowed
    assert q.overlaps_any(s, maxgap=0).tolist() == [True, True, True, True]


def test_many_to_one_large() -> None:
    # exercises chunking path through the candidate expansion
    rng = np.random.default_rng(0)
    qs = rng.integers(1, 1000, 500)
    x = IRanges(start=qs, width=rng.integers(1, 50, 500))
    h = x.find_overlaps(x)
    assert h.count_per_query().min() >= 1
    assert len(h) == int(x.count_overlaps(x).sum())


def brute(q: IRanges, s: IRanges, type: str, maxgap: int, minoverlap: int) -> list[tuple[int, int]]:
    """Independent per-pair oracle (python loops)."""
    out = []
    for i in range(len(q)):
        qs, qe = int(q.start[i]), int(q.end[i])
        for j in range(len(s)):
            ss, se = int(s.start[j]), int(s.end[j])
            m = max(maxgap, 0)
            if type == "any":
                if qe < qs and se < ss:
                    ok = max(qs, ss) - min(qe, se) - 1 <= maxgap
                elif qe < qs:  # empty query at qs
                    ok = (ss < qs <= se) or (max(qs, ss) - min(qe, se) - 1 <= maxgap)
                elif se < ss:
                    ok = (qs < ss <= qe) or (max(qs, ss) - min(qe, se) - 1 <= maxgap)
                else:
                    ok = max(qs, ss) <= min(qe, se) + 1 + maxgap
            elif type == "start":
                ok = abs(qs - ss) <= m
            elif type == "end":
                ok = abs(qe - se) <= m
            elif type == "equal":
                ok = abs(qs - ss) <= m and abs(qe - se) <= m
            else:
                ok = qs >= ss and qe <= se and (maxgap < 0 or (qs - ss) + (se - qe) <= maxgap)
            if ok and minoverlap > 0:
                ok = min(qe, se) - max(qs, ss) + 1 >= minoverlap
            if ok:
                out.append((i, j))
    return out


@given(
    iranges_st(max_size=10),
    iranges_st(max_size=10),
    st.sampled_from(["any", "start", "end", "within", "equal"]),
    st.integers(-1, 6),
    st.integers(0, 4),
)
def test_matches_brute_force(q: IRanges, s: IRanges, type: str, maxgap: int, minoverlap: int) -> None:
    if type == "any" and maxgap != -1 and minoverlap != 0:
        minoverlap = 0
    h = q.find_overlaps(s, maxgap=maxgap, minoverlap=minoverlap, type=type)
    assert h.pairs() == brute(q, s, type, maxgap, minoverlap)
    first = q.find_overlaps(s, maxgap=maxgap, minoverlap=minoverlap, type=type, select="first")
    expect = {}
    for i, j in h.pairs():
        expect.setdefault(i, j)
    assert first.tolist() == [expect.get(i, NO_HIT) for i in range(len(q))]
