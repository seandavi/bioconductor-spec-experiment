"""GRanges / Seqinfo unit tests beyond the shared cases."""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e4 import GRanges, Seqinfo


def test_seqinfo_merge_is_symmetric_on_values() -> None:
    a = Seqinfo(["x", "y"], [10, None], [None, True], "g")
    b = Seqinfo(["y", "z"], [20, 30], [None, None], None)
    ab, ba = a.merge(b), b.merge(a)
    assert ab.seqlengths == ba.seqlengths | {} and ab.seqnames == ("x", "y", "z")
    assert ba.seqnames == ("y", "z", "x")
    assert ab.genome == {"x": "g", "y": "g", "z": "g"} or ab.genome["x"] == "g"


def test_seqinfo_is_immutable() -> None:
    si = Seqinfo(["a"])
    with pytest.raises(AttributeError):
        si.seqnames = ("b",)  # type: ignore[misc]


def test_granges_getitem_does_not_alias_input() -> None:
    g = GRanges(["a", "b"], start=[1, 2], width=[1, 1])
    sub = g[[0]]
    assert len(sub) == 1 and len(g) == 2
    assert sub.seqinfo is g.seqinfo


def test_trim_leaves_in_bounds_ranges_alone() -> None:
    g = GRanges(["a"], start=[2], width=[3], seqinfo=Seqinfo(["a"], [10]))
    t = g.trim()
    assert t.start.tolist() == [2] and t.width.tolist() == [3]


def test_concat_remaps_codes_by_name() -> None:
    a = GRanges(["x"], start=[1], width=[1], seqinfo=Seqinfo(["x", "y"]))
    b = GRanges(["y"], start=[1], width=[1], seqinfo=Seqinfo(["y", "x"]))
    c = GRanges.concat(a, b)
    assert c.seqnames.tolist() == ["x", "y"] and c.seqlevels == ["x", "y"]


def test_strand_codes_roundtrip() -> None:
    g = GRanges(["a", "a", "a"], start=[1, 2, 3], width=[1, 1, 1], strand=["+", "-", "*"])
    assert g.strand.tolist() == ["+", "-", "*"]
    assert np.array_equal(g._strand, [0, 1, 2])


def _random_granges(seed: int, n: int = 14) -> GRanges:
    rng = np.random.default_rng(seed)
    return GRanges(
        rng.choice(["c1", "c2", "c3"], n).tolist(),
        start=rng.integers(1, 40, n),
        width=rng.integers(0, 8, n),
        strand=rng.choice(["+", "-", "*"], n).tolist(),
        seqinfo=Seqinfo(["c1", "c2", "c3"], [60, 60, 60]),
    )


def _cover(g: GRanges) -> set[tuple[str, str, int]]:
    return {
        (sq, st, p)
        for sq, st, s, w in zip(g.seqnames, g.strand, g.start.tolist(), g.width.tolist(), strict=True)
        for p in range(s, s + w)
    }


@pytest.mark.parametrize("seed", range(15))
def test_group_ops_preserve_per_group_coverage(seed: int) -> None:
    g = _random_granges(seed)
    assert _cover(g.reduce()) == _cover(g)
    assert _cover(g.disjoin()) == _cover(g)
    # gaps + coverage tile every (seqlevel, strand) group exactly
    full = {(sq, st, p) for sq in ("c1", "c2", "c3") for st in "+-*" for p in range(1, 61)}
    assert _cover(g.gaps()) == full - _cover(g)


@pytest.mark.parametrize("seed", range(15))
def test_ignore_strand_pools_strands(seed: int) -> None:
    g = _random_granges(seed)
    r = g.reduce(ignore_strand=True)
    assert set(r.strand.tolist()) <= {"*"}
    pooled = {(sq, p) for sq, _, p in _cover(g)}
    assert {(sq, p) for sq, _, p in _cover(r)} == pooled


@pytest.mark.parametrize("seed", range(15))
def test_set_ops_match_per_group_set_algebra(seed: int) -> None:
    x, y = _random_granges(seed), _random_granges(seed + 99)
    cx, cy = _cover(x), _cover(y)
    assert _cover(x.union(y)) == cx | cy
    assert _cover(x.intersect(y)) == cx & cy
    assert _cover(x.setdiff(y)) == cx - cy
