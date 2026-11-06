"""Invariants that need no reference implementation (hypothesis)."""

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from ranges_e0 import GRanges, IRanges

SETTINGS = settings(max_examples=150, deadline=None)

iv = st.tuples(st.integers(-3, 30), st.integers(0, 8))
irange_lists = st.lists(iv, min_size=0, max_size=12)


def make(spec) -> IRanges:
    return IRanges([a for a, _ in spec], width=[w for _, w in spec]) if spec else IRanges.empty()


def positions(x: IRanges) -> set[int]:
    return {p for s, e in x.to_tuples() for p in range(s, e + 1)}


def gen_positions(x: IRanges) -> list[int]:
    return sorted(positions(x))


@SETTINGS
@given(irange_lists, st.integers(0, 4), st.booleans())
def test_reduce_invariants(spec, mgw, drop):
    x = make(spec)
    r = x.reduce(min_gapwidth=mgw, drop_empty_ranges=drop)
    t = r.to_tuples()
    # sorted by start, and consecutive ranges are separated by at least min_gapwidth
    for (s1, e1), (s2, _) in zip(t, t[1:]):
        assert s2 - e1 - 1 >= mgw
    if mgw <= 1:
        assert positions(r) == positions(x)  # same covered positions
    else:
        # a wider min.gapwidth bridges gaps (and nearby empty ranges merge into
        # non-empty ones), so the result may cover more than the input
        assert positions(r) >= positions(x)
    again = r.reduce(min_gapwidth=mgw, drop_empty_ranges=drop)
    assert again == r  # idempotent
    if drop:
        assert all(w > 0 for w in r.width)


@SETTINGS
@given(irange_lists)
def test_reduce_revmap_partitions_the_input(spec):
    x = make(spec)
    r = x.reduce(with_revmap=True)
    members = [i for grp in r.mcols["revmap"].to_list() for i in grp]
    assert sorted(members) == list(range(len(x)))
    for (s, e), grp in zip(r.to_tuples(), r.mcols["revmap"].to_list()):
        assert s == min(x.start[grp]) and e == max(x.end[grp])


@SETTINGS
@given(irange_lists, irange_lists)
def test_set_ops_match_position_algebra(a, b):
    x, y = make(a), make(b)
    px, py = positions(x), positions(y)
    assert positions(x.union(y)) == px | py
    assert positions(x.intersect(y)) == px & py
    assert positions(x.setdiff(y)) == px - py
    for r in (x.union(y), x.intersect(y), x.setdiff(y)):
        assert r == r.reduce(drop_empty_ranges=True)  # canonical form


@SETTINGS
@given(irange_lists)
def test_disjoin_and_gaps(spec):
    x = make(spec)
    d = x.disjoin()
    assert positions(d) == positions(x)
    assert int(d.width.sum()) == len(positions(x))  # pieces never overlap
    r = x.reduce(drop_empty_ranges=True)
    if len(r):
        lo, hi = int(r.start[0]), int(r.end[-1])
        g = x.gaps()
        assert positions(g) | positions(r) == set(range(lo, hi + 1))
        assert not positions(g) & positions(r)
        g2 = x.gaps(start=lo - 2, end=hi + 3)
        assert positions(g2) | positions(r) == set(range(lo - 2, hi + 4)) | positions(r)
        assert not positions(g2) & positions(r)


@SETTINGS
@given(irange_lists)
def test_sort_unique_rank(spec):
    x = make(spec)
    assert x.sort() == x[x.order()]
    keys = list(zip(x.sort().start.tolist(), x.sort().width.tolist()))
    assert keys == sorted(keys)
    assert sorted(x.rank("first").tolist()) == list(range(1, len(x) + 1))
    u = x.unique()
    assert len(u) == len(set(zip(x.start.tolist(), x.width.tolist())))
    assert not u.duplicated().any()


@SETTINGS
@given(irange_lists, st.integers(-4, 6), st.integers(0, 10))
def test_intra_range_properties(spec, sh, w):
    x = make(spec)
    assert x.shift(sh).shift(-sh) == x
    assert x.resize(w).width.tolist() == [w] * len(x)
    assert x.resize(w, fix="end").end.tolist() == x.end.tolist()
    assert x.resize(w).start.tolist() == x.start.tolist()
    assert x.flank(w).end.tolist() == (x.start - 1).tolist()
    assert x.flank(w, start=False).start.tolist() == (x.end + 1).tolist()
    assert x.flank(w, both=True).width.tolist() == [2 * w] * len(x)


@SETTINGS
@given(irange_lists)
def test_distance_is_symmetric_and_zero_when_overlapping(spec):
    x = make(spec)
    y = x[::-1]
    assert (x.distance(y) == y.distance(x)[::-1]).all()
    # R (2.47.3): empty ranges never overlap, so an empty range is at distance 1 from itself
    assert (x.distance(x) == (x.width == 0).astype(int)).all()


@SETTINGS
@given(irange_lists)
def test_distance_is_gap_plus_one_unless_overlapping(spec):
    x = make(spec)
    y = x[::-1]
    lo = np.maximum(x.start, y.start)
    hi = np.minimum(x.end, y.end)
    assert (x.distance(y) == np.where(lo <= hi, 0, lo - hi)).all()


# ---- GRanges ----------------------------------------------------------------------
g_row = st.tuples(st.sampled_from(["chr1", "chr2", "chr3"]), st.integers(1, 30), st.integers(0, 6),
                  st.sampled_from(["+", "-", "*"]))
g_lists = st.lists(g_row, min_size=0, max_size=14)


def gmake(spec) -> GRanges:
    if not spec:
        return GRanges.empty()
    return GRanges([c for c, *_ in spec], IRanges([s for _, s, _, _ in spec], width=[w for _, _, w, _ in spec]),
                   [t for *_, t in spec])


def gpositions(g: GRanges, ignore_strand: bool = False) -> set[tuple[str, str, int]]:
    out = set()
    for c, s, e, t in g.to_tuples():
        for p in range(s, e + 1):
            out.add((c, "*" if ignore_strand else t, p))
    return out


@SETTINGS
@given(g_lists, g_lists, st.booleans())
def test_granges_set_ops_match_position_algebra(a, b, ignore):
    x, y = gmake(a), gmake(b)
    px, py = gpositions(x, ignore), gpositions(y, ignore)
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert gpositions(x.union(y, ignore_strand=ignore)) == px | py
        assert gpositions(x.intersect(y, ignore_strand=ignore)) == px & py
        assert gpositions(x.setdiff(y, ignore_strand=ignore)) == px - py


@SETTINGS
@given(g_lists, st.booleans())
def test_granges_reduce_and_disjoin_keep_positions(a, ignore):
    x = gmake(a)
    assert gpositions(x.reduce(ignore_strand=ignore)) == gpositions(x, ignore)
    assert gpositions(x.disjoin(ignore_strand=ignore)) == gpositions(x, ignore)
    r = x.reduce(ignore_strand=ignore)
    assert r.order().tolist() == list(range(len(r)))  # already in (seqlevel, strand, start) order


@SETTINGS
@given(g_lists, g_lists, st.booleans(), st.integers(-1, 3))
def test_granges_overlaps_match_bruteforce(a, b, ignore, maxgap):
    q, s = gmake(a), gmake(b)
    got = [] if not len(q) or not len(s) else list(q.find_overlaps(s, maxgap=maxgap, ignore_strand=ignore))
    exp = []
    for i, (qc, qs, qe, qt) in enumerate(q.to_tuples() if len(q) else []):
        for j, (sc, ss, se, stt) in enumerate(s.to_tuples() if len(s) else []):
            if qc != sc:
                continue
            if not ignore and not (qt == "*" or stt == "*" or qt == stt):
                continue
            if qs <= se + maxgap + 1 and ss <= qe + maxgap + 1:
                exp.append((i, j))
    assert got == exp


@SETTINGS
@given(g_lists, st.booleans())
def test_granges_sort_is_sorted_by_code_keys(a, ignore):
    x = gmake(a)
    if not len(x):
        return
    s = x.sort(ignore_strand=ignore)
    lv = {n: i for i, n in enumerate(x.seqlevels)}
    sv = {"+": 0, "-": 1, "*": 2}
    keys = [(lv[c], 0 if ignore else sv[t], st_, e - st_ + 1) for c, st_, e, t in s.to_tuples()]
    assert keys == sorted(keys)


@SETTINGS
@given(g_lists)
def test_granges_coverage_total_is_sum_of_widths(a):
    x = gmake(a)
    cov = x.coverage()
    total = sum(int((r.values * r.lengths).sum()) for r in cov.values())
    assert total == int(x.width.sum())
