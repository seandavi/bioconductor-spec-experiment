"""Randomised checks of the array kernels against brute-force definitions of the same rules."""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e2 import _inter, _overlaps
from ranges_e2._nearest import nearest_ties, pair_distance

RNG = np.random.default_rng(20261105)


def _random_ranges(n: int, span: int = 30, max_width: int = 6) -> tuple[np.ndarray, np.ndarray]:
    return RNG.integers(-3, span, n), RNG.integers(0, max_width + 1, n)


def _positions(s: int, w: int) -> set[int]:
    return set(range(s, s + w))


@pytest.mark.parametrize("seed", range(40))
def test_reduce_matches_position_union(seed: int) -> None:
    s, w = _random_ranges(int(RNG.integers(0, 9)))
    rs, rw, _ = _inter.reduce_ranges(s, w, 1, drop_empty=True)
    covered: set[int] = set().union(*[_positions(a, b) for a, b in zip(s, w)])
    got: set[int] = set().union(*[_positions(a, b) for a, b in zip(rs, rw)]) if len(rs) else set()
    assert got == covered
    assert all(rs[i] + rw[i] < rs[i + 1] for i in range(len(rs) - 1))  # separated by >= 1 position


@pytest.mark.parametrize("seed", range(40))
def test_gaps_and_setops_match_position_sets(seed: int) -> None:
    xs, xw = _random_ranges(int(RNG.integers(0, 7)))
    ys, yw = _random_ranges(int(RNG.integers(0, 7)))
    cx = set().union(*[_positions(a, b) for a, b in zip(xs, xw)])
    cy = set().union(*[_positions(a, b) for a, b in zip(ys, yw)])

    def cover(s: np.ndarray, w: np.ndarray) -> set[int]:
        return set().union(*[_positions(a, b) for a, b in zip(s, w)]) if len(s) else set()

    assert cover(*_inter.union_ranges(xs, xw, ys, yw)) == cx | cy
    assert cover(*_inter.intersect_ranges(xs, xw, ys, yw)) == cx & cy
    assert cover(*_inter.setdiff_ranges(xs, xw, ys, yw)) == cx - cy
    lo, hi = -3, 33
    gs, gw = _inter.gaps_ranges(xs, xw, lo, hi)
    assert cover(gs, gw) == set(range(lo, hi + 1)) - cx


@pytest.mark.parametrize("seed", range(30))
@pytest.mark.parametrize("kind", ["any", "start", "end", "within", "equal"])
def test_overlap_pairs_match_brute_force(seed: int, kind: str) -> None:
    qs, qw = _random_ranges(int(RNG.integers(0, 8)))
    ss, sw = _random_ranges(int(RNG.integers(0, 8)))
    maxgap = int(RNG.choice([-1, 0, 1, 3])) if kind != "any" else -1
    minoverlap = int(RNG.choice([0, 1, 2]))
    q, s = _overlaps.overlap_pairs(qs, qw, ss, sw, maxgap, minoverlap, kind)
    got = set(zip(q.tolist(), s.tolist()))
    want = set()
    m = max(maxgap, 0)
    for i in range(len(qs)):
        for j in range(len(ss)):
            qe, se = qs[i] + qw[i] - 1, ss[j] + sw[j] - 1
            shared = min(qe, se) - max(qs[i], ss[j]) + 1
            if kind == "any":
                ok = qs[i] <= se and ss[j] <= qe
            elif kind == "start":
                ok = abs(qs[i] - ss[j]) <= m
            elif kind == "end":
                ok = abs(qe - se) <= m
            elif kind == "equal":
                ok = abs(qs[i] - ss[j]) <= m and abs(qe - se) <= m
            else:
                ok = qs[i] >= ss[j] and qe <= se and (maxgap <= 0 or sw[j] - qw[i] <= maxgap)
            if ok and (minoverlap == 0 or shared >= minoverlap):
                want.add((i, j))
    assert got == want
    assert q.tolist() == sorted(q.tolist())


@pytest.mark.parametrize("seed", range(30))
def test_nearest_ties_are_minimum_distance(seed: int) -> None:
    qs, qw = _random_ranges(int(RNG.integers(1, 6)))
    ss, sw = _random_ranges(int(RNG.integers(1, 8)))
    for i, ties in enumerate(nearest_ties(qs, qw, ss, sw)):
        d = pair_distance(qs[i], qs[i] + qw[i] - 1, ss, ss + sw - 1)
        assert set(ties.tolist()) == set(np.nonzero(d == d.min())[0].tolist())


def test_distance_adjacent_is_one_since_iranges_2_47_3() -> None:
    from ranges_e2 import IRanges

    x = IRanges(start=[1, 1, 1, 1], width=[5, 5, 5, 5])
    y = IRanges(start=[3, 6, 7, 9], width=[6, 3, 2, 4])
    assert x.distance(y).tolist() == [0, 1, 2, 4]
    assert y.distance(x).tolist() == [0, 1, 2, 4]
