import warnings

import pytest
from hypothesis import given, strategies as st

from helpers import se
from ranges_e0b import GRanges, IRanges, Seqinfo


def rows(x: GRanges) -> list[tuple[str, int, int, str]]:
    return list(zip(x.seqnames, x.start.tolist(), x.end.tolist(), x.strand))


def sample() -> GRanges:
    return GRanges(
        ["chr1", "chr1", "chr1", "chr2", "chr1"],
        IRanges(start=[1, 5, 5, 3, 20], end=[10, 12, 8, 6, 25]),
        strand=["+", "+", "-", "+", "+"],
    )


def test_reduce_per_seq_and_strand() -> None:
    # R: reduce(gr)
    assert rows(sample().reduce()) == [
        ("chr1", 1, 12, "+"), ("chr1", 20, 25, "+"), ("chr1", 5, 8, "-"), ("chr2", 3, 6, "+")
    ]


def test_reduce_ignore_strand() -> None:
    assert rows(sample().reduce(ignore_strand=True)) == [
        ("chr1", 1, 12, "*"), ("chr1", 20, 25, "*"), ("chr2", 3, 6, "*")
    ]


def test_reduce_revmap_and_drop_names() -> None:
    x = sample().with_names(list("abcde")).with_mcols({"v": [1, 2, 3, 4, 5]})
    r, rev = x.reduce(with_revmap=True)
    assert [m.tolist() for m in rev] == [[0, 1], [4], [2], [3]]
    assert r.names is None and r.mcols.width == 0


def test_reduce_min_gapwidth() -> None:
    x = GRanges("a", IRanges(start=[1, 12], end=[10, 20]))
    assert len(x.reduce()) == 2
    assert len(x.reduce(min_gapwidth=2)) == 1


def test_range() -> None:
    assert rows(sample().range()) == [
        ("chr1", 1, 25, "+"), ("chr1", 5, 8, "-"), ("chr2", 3, 6, "+")
    ]
    assert rows(sample().range(ignore_strand=True)) == [("chr1", 1, 25, "*"), ("chr2", 3, 6, "*")]


def test_disjoin() -> None:
    assert rows(sample().disjoin()) == [
        ("chr1", 1, 4, "+"), ("chr1", 5, 10, "+"), ("chr1", 11, 12, "+"), ("chr1", 20, 25, "+"),
        ("chr1", 5, 8, "-"), ("chr2", 3, 6, "+"),
    ]
    assert rows(sample().disjoin(ignore_strand=True))[:4] == [
        ("chr1", 1, 4, "*"), ("chr1", 5, 8, "*"), ("chr1", 9, 10, "*"), ("chr1", 11, 12, "*"),
    ]


def test_gaps_all_three_strands_with_seqlengths() -> None:
    # R: gaps(GRanges("chr1", IRanges(5,10), "+", seqinfo=Seqinfo("chr1", 100)))
    x = GRanges("chr1", IRanges(start=5, end=10), "+", seqinfo=Seqinfo(["chr1"], [100]))
    assert rows(x.gaps()) == [
        ("chr1", 1, 4, "+"), ("chr1", 11, 100, "+"), ("chr1", 1, 100, "-"), ("chr1", 1, 100, "*")
    ]
    assert rows(x.gaps(ignore_strand=True)) == [("chr1", 1, 4, "*"), ("chr1", 11, 100, "*")]


def test_gaps_unknown_seqlength_uses_max_end_and_skips_empty_levels() -> None:
    # UNVERIFIED: R behaviour with NA seqlengths
    x = GRanges("chr1", IRanges(start=5, end=10), "+", seqinfo=Seqinfo(["chr1", "chr2"]))
    assert rows(x.gaps()) == [
        ("chr1", 1, 4, "+"), ("chr1", 1, 10, "-"), ("chr1", 1, 10, "*")
    ]


def test_gaps_empty_seqlevel_with_length() -> None:
    x = GRanges("chr1", IRanges(start=5, end=10), "+", seqinfo=Seqinfo(["chr1", "chr2"], [20, 30]))
    out = rows(x.gaps(ignore_strand=True))
    assert out == [("chr1", 1, 4, "*"), ("chr1", 11, 20, "*"), ("chr2", 1, 30, "*")]


def test_gaps_explicit_bounds() -> None:
    x = GRanges("chr1", IRanges(start=5, end=10), "*")
    assert rows(x.gaps(start=3, end=12, ignore_strand=True)) == [("chr1", 3, 4, "*"), ("chr1", 11, 12, "*")]
    assert rows(x.gaps(start={"chr1": 7}, end={"chr1": 20}, ignore_strand=True)) == [("chr1", 11, 20, "*")]


def test_setops_strand_aware() -> None:
    a = GRanges(["c", "c"], IRanges(start=[1, 1], end=[10, 10]), ["+", "-"])
    b = GRanges(["c", "c"], IRanges(start=[5, 5], end=[20, 7]), ["+", "+"])
    assert rows(a.union(b)) == [("c", 1, 20, "+"), ("c", 1, 10, "-")]
    assert rows(a.intersect(b)) == [("c", 5, 10, "+")]
    assert rows(a.setdiff(b)) == [("c", 1, 4, "+"), ("c", 1, 10, "-")]
    assert rows(a.union(b, ignore_strand=True)) == [("c", 1, 20, "*")]
    assert rows(a.intersect(b, ignore_strand=True)) == [("c", 5, 10, "*")]
    assert rows(a.setdiff(b, ignore_strand=True)) == [("c", 1, 4, "*")]


def test_setops_merge_seqinfo_and_levels() -> None:
    a = GRanges("x", IRanges(start=1, end=5), "+", seqinfo=Seqinfo(["x"], [50]))
    b = GRanges("y", IRanges(start=1, end=5), "+", seqinfo=Seqinfo(["y"], [60]))
    with pytest.warns(UserWarning, match="no sequence names in common"):
        u = a.union(b)
    assert rows(u) == [("x", 1, 5, "+"), ("y", 1, 5, "+")]
    assert u.seqlengths == {"x": 50, "y": 60}


def test_parallel_ops() -> None:
    a = GRanges(["c", "c"], IRanges(start=[1, 1], end=[10, 10], names=["p", "q"]), ["+", "*"], mcols={"v": [1, 2]})
    b = GRanges(["c", "c"], IRanges(start=[5, 5], end=[20, 20]), ["*", "-"])
    p = a.pintersect(b)
    assert rows(p) == [("c", 5, 10, "+"), ("c", 5, 10, "-")] and p.names == ["p", "q"]
    assert p.mcols["v"].to_list() == [1, 2]
    assert rows(a.punion(b)) == [("c", 1, 20, "+"), ("c", 1, 20, "-")]
    assert rows(a.pintersect(b, ignore_strand=True))[0][3] == "*"
    assert rows(a.psetdiff(b)) == [("c", 1, 4, "+"), ("c", 1, 4, "*")]


def test_parallel_errors() -> None:
    a = GRanges("c", IRanges(start=1, end=5), "+")
    with pytest.raises(ValueError, match="strands"):
        a.pintersect(GRanges("c", IRanges(start=1, end=5), "-"))
    with pytest.warns(UserWarning, match="no sequence names"), pytest.raises(ValueError, match="seqnames"):
        a.pintersect(GRanges("d", IRanges(start=1, end=5), "+"))
    with pytest.raises(ValueError, match="lengths"):
        a.punion(a.concat(a))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert len(a.pintersect(GRanges("c", IRanges(start=1, end=5), "-"), ignore_strand=True)) == 1


@given(
    st.lists(st.tuples(st.sampled_from(["a", "b"]), st.sampled_from("+-*"), st.integers(0, 20), st.integers(0, 6)), max_size=15),
    st.lists(st.tuples(st.sampled_from(["a", "b"]), st.sampled_from("+-*"), st.integers(0, 20), st.integers(0, 6)), max_size=15),
    st.booleans(),
)
def test_granges_ops_equal_iranges_per_group(xs: list, ys: list, ignore: bool) -> None:  # type: ignore[type-arg]
    info = Seqinfo(["a", "b"])

    def mk(rs: list) -> GRanges:  # type: ignore[type-arg]
        return GRanges([r[0] for r in rs], IRanges(start=[r[2] for r in rs], width=[r[3] for r in rs]) if rs else IRanges(),
                       [r[1] for r in rs], seqinfo=info)

    x, y = mk(xs), mk(ys)

    def sub(rs: list, key: tuple[str, str]) -> IRanges:  # type: ignore[type-arg]
        sel = [r for r in rs if r[0] == key[0] and (ignore or r[1] == key[1])]
        return IRanges(start=[r[2] for r in sel], width=[r[3] for r in sel]) if sel else IRanges()

    keys = [("a", "+"), ("a", "-"), ("a", "*"), ("b", "+"), ("b", "-"), ("b", "*")]
    if ignore:
        keys = [("a", "*"), ("b", "*")]
    for op in ("reduce", "disjoin", "range"):
        got = getattr(x, op)(ignore_strand=ignore)
        expect = []
        for key in keys:
            part = getattr(sub(xs, key), op)()
            expect += [(key[0], s, e, key[1]) for s, e in se(part)]
        assert rows(got) == expect, op
    for op in ("union", "intersect", "setdiff"):
        got = getattr(x, op)(y, ignore_strand=ignore)
        expect = []
        for key in keys:
            part = getattr(sub(xs, key), op)(sub(ys, key))
            expect += [(key[0], s, e, key[1]) for s, e in se(part)]
        assert rows(got) == expect, op
