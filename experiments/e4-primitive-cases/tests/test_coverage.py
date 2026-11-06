"""Coverage vs a naive per-position count on seeded random inputs."""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e4 import IRanges, Rle


@pytest.mark.parametrize("seed", range(25))
def test_coverage_matches_naive_count(seed: int) -> None:
    rng = np.random.default_rng(seed)
    x = IRanges(start=rng.integers(-3, 30, 10), width=rng.integers(0, 8, 10))
    w = rng.integers(1, 4, 10)
    cov = x.coverage(weight=w).decode()
    length = max(int(x.end.max()), 0)
    naive = np.zeros(length, dtype=np.int64)
    for s, wd, wt in zip(x.start.tolist(), x.width.tolist(), w.tolist(), strict=True):
        for p in range(s, s + wd):
            if 1 <= p <= length:
                naive[p - 1] += wt
    assert cov.tolist() == naive.tolist()


def test_rle_encode_decode_roundtrip() -> None:
    v = np.array([0, 0, 1, 1, 1, 0, 2])
    r = Rle.encode(v)
    assert r.values.tolist() == [0, 1, 0, 2] and r.lengths.tolist() == [2, 3, 1, 1]
    assert r.decode().tolist() == v.tolist() and len(r) == 7


def test_rle_encode_empty() -> None:
    r = Rle.encode(np.zeros(0, dtype=np.int64))
    assert len(r) == 0 and r.values.size == 0
