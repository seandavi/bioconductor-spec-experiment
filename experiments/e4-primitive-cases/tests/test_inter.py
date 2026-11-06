"""Brute-force oracles for inter-range kernels on seeded random inputs."""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e4 import IRanges

rng_seeds = range(25)


def random_ranges(seed: int, n: int = 12, span: int = 40, max_w: int = 8) -> IRanges:
    rng = np.random.default_rng(seed)
    return IRanges(start=rng.integers(-5, span, n), width=rng.integers(0, max_w, n))


def coverage_set(x: IRanges) -> set[int]:
    return {p for s, w in zip(x.start.tolist(), x.width.tolist(), strict=True) for p in range(s, s + w)}


@pytest.mark.parametrize("seed", rng_seeds)
def test_reduce_default_covers_same_positions(seed: int) -> None:
    x = random_ranges(seed)
    r = x.reduce()
    assert coverage_set(r) == coverage_set(x)
    nonempty = r[r.width > 0]
    s, e = nonempty.start, nonempty.end
    assert (s[1:] > e[:-1] + 1).all()  # disjoint and non-adjacent


@pytest.mark.parametrize("seed", rng_seeds)
def test_reduce_revmap_partitions_input(seed: int) -> None:
    x = random_ranges(seed)
    r = x.reduce(with_revmap=True, min_gapwidth=2)
    members = [m for row in r.mcols["revmap"].to_list() for m in row]
    assert sorted(members) == list(range(len(x)))


@pytest.mark.parametrize("seed", rng_seeds)
def test_disjoin_pieces_have_constant_cover(seed: int) -> None:
    x = random_ranges(seed)
    d = x.disjoin(with_revmap=True)
    assert coverage_set(d) == coverage_set(x)
    for s, w, who in zip(d.start.tolist(), d.width.tolist(), d.mcols["revmap"].to_list(), strict=True):
        for p in range(s, s + w):
            covering = [
                i for i, (xs, xw) in enumerate(zip(x.start.tolist(), x.width.tolist(), strict=True)) if xs <= p < xs + xw
            ]
            assert covering == who


@pytest.mark.parametrize("seed", rng_seeds)
def test_gaps_is_complement_within_bounds(seed: int) -> None:
    x = random_ranges(seed)
    lo, hi = -5, 50
    g = x.gaps(start=lo, end=hi)
    assert coverage_set(g) == set(range(lo, hi + 1)) - coverage_set(x)
    assert (g.width > 0).all()


def test_reduce_min_gapwidth_rejects_negative() -> None:
    with pytest.raises(ValueError):
        IRanges(start=[1], width=[1]).reduce(min_gapwidth=-1)
