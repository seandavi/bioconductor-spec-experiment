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


def _two_seq_granges(circular: bool = False) -> "GRanges":
    from ranges_e4 import GRanges, Seqinfo

    si = Seqinfo(["a", "b"], [10, 10], [circular, circular])
    return GRanges(["a", "b"], start=[1, 2], width=[3, 3], seqinfo=si)


@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize("weight", [[1, 2, 3], [1, 2, 3, 4], []])
def test_granges_coverage_weight_length_mismatch_raises(circular: bool, weight: list[int]) -> None:
    with pytest.raises(ValueError):
        _two_seq_granges(circular).coverage(weight=weight)


@pytest.mark.parametrize("circular", [False, True])
def test_granges_coverage_valid_weight_vector_is_split_per_seqlevel(circular: bool) -> None:
    cov = _two_seq_granges(circular).coverage(weight=[5, 7])
    assert cov["a"].decode().tolist() == [5, 5, 5] + [0] * 7
    assert cov["b"].decode().tolist() == [0, 7, 7, 7] + [0] * 6


def test_iranges_coverage_weight_length_mismatch_raises() -> None:
    x = IRanges(start=[1, 2], width=[1, 1])
    for bad in ([1, 2, 3], []):
        with pytest.raises(ValueError):
            x.coverage(weight=bad)
