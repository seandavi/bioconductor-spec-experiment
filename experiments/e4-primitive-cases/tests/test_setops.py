"""Set-op oracles: compare against Python sets of covered positions."""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e4 import IRanges


def positions(x: IRanges) -> set[int]:
    return {p for s, w in zip(x.start.tolist(), x.width.tolist(), strict=True) for p in range(s, s + w)}


def rand(seed: int, n: int = 8) -> IRanges:
    rng = np.random.default_rng(seed)
    return IRanges(start=rng.integers(0, 30, n), width=rng.integers(0, 7, n))


@pytest.mark.parametrize("seed", range(30))
def test_sequence_ops_match_set_algebra(seed: int) -> None:
    x, y = rand(seed), rand(seed + 1000)
    px, py = positions(x), positions(y)
    for result, expected in (
        (x.union(y), px | py),
        (x.intersect(y), px & py),
        (x.setdiff(y), px - py),
    ):
        assert positions(result) == expected
        assert (result.width > 0).all()
        # normalised: sorted and non-adjacent
        assert (result.start[1:] > result.end[:-1] + 1).all()


@pytest.mark.parametrize("seed", range(30))
def test_pintersect_pairwise(seed: int) -> None:
    x, y = rand(seed), rand(seed + 500)
    r = x.pintersect(y, resolve_empty="start")
    for i in range(len(x)):
        px, py = positions(x[[i]]), positions(y[[i]])
        assert positions(r[[i]]) == px & py


def test_parallel_ops_reject_unequal_lengths() -> None:
    a, b = IRanges(start=[1, 2], width=[1, 1]), IRanges(start=[1], width=[1])
    for f in (a.punion, a.pintersect, a.psetdiff):
        with pytest.raises(ValueError):
            f(b)
