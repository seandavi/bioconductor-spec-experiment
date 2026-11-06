"""Ordering invariants on seeded random inputs."""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e4 import IRanges


def rand(seed: int) -> IRanges:
    rng = np.random.default_rng(seed)
    return IRanges(start=rng.integers(0, 6, 15), width=rng.integers(0, 4, 15))


@pytest.mark.parametrize("seed", range(20))
def test_order_matches_python_sorted(seed: int) -> None:
    x = rand(seed)
    keys = list(zip(x.start.tolist(), x.width.tolist(), strict=True))
    assert x.order().tolist() == sorted(range(len(x)), key=lambda i: keys[i])
    assert x.order(decreasing=True).tolist() == sorted(range(len(x)), key=lambda i: (-keys[i][0], -keys[i][1]))


@pytest.mark.parametrize("seed", range(20))
def test_rank_methods_consistent(seed: int) -> None:
    x = rand(seed)
    mn, mx, av = (x.rank(ties_method=m) for m in ("min", "max", "average"))
    assert ((mn + mx) / 2 == av).all()
    first, last = x.rank(ties_method="first"), x.rank(ties_method="last")
    assert sorted(first.tolist()) == list(range(1, len(x) + 1))
    assert sorted(last.tolist()) == list(range(1, len(x) + 1))
    assert (mn <= first).all() and (first <= mx).all()


@pytest.mark.parametrize("seed", range(20))
def test_duplicated_and_unique_agree(seed: int) -> None:
    x = rand(seed)
    keys = list(zip(x.start.tolist(), x.width.tolist(), strict=True))
    seen: set[tuple[int, int]] = set()
    expect = []
    for k in keys:
        expect.append(k in seen)
        seen.add(k)
    assert x.duplicated().tolist() == expect
    assert len(x.unique()) == len(set(keys))
    last_expect = [k in set(keys[i + 1 :]) for i, k in enumerate(keys)]
    assert x.duplicated(from_last=True).tolist() == last_expect


def test_rank_requires_ties_method() -> None:
    with pytest.raises(TypeError):
        rand(0).rank()  # type: ignore[call-arg]
