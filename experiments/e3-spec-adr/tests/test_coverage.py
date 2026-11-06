"""coverage and Rle (spec/50-coverage.md)."""

from __future__ import annotations

import numpy as np
import pytest

from helpers import ir
from ranges_e3.rle import Rle


def runs(r: Rle) -> list[tuple[int, int]]:
    return list(zip(r.lengths.tolist(), r.values.tolist()))


@pytest.mark.spec("COV.1", "COV.2")
def test_coverage_basic_and_zero_width() -> None:
    r = ir((1, 3), (5, 4), (10, 13), (10, 10)).coverage()
    assert runs(r) == [(3, 1), (6, 0), (1, 2), (3, 1)]
    assert runs(ir((5, 4)).coverage()) == [(4, 0)]
    assert runs(ir().coverage()) == []


@pytest.mark.spec("COV.1", "COV.2")
def test_coverage_width_shift_negative() -> None:
    x = ir((2, 3), (4, 6))
    assert runs(x.coverage(width=10)) == [(1, 0), (5, 1), (4, 0)]
    assert runs(x.coverage(width=4)) == [(1, 0), (3, 1)]
    assert runs(x.coverage(shift=1)) == [(2, 0), (5, 1)]
    assert runs(x.coverage(shift=-1)) == [(5, 1)]
    assert runs(x.coverage(shift=2, width=6)) == [(3, 0), (3, 1)]
    assert runs(ir((-2, 3)).coverage()) == [(3, 1)]
    assert runs(ir((-5, -3)).coverage()) == []
    assert runs(ir((1, 2), (1, 2)).coverage(shift=[0, 5])) == [(2, 1), (3, 0), (2, 1)]
    assert runs(ir().coverage(width=5)) == [(5, 0)]
    assert runs(ir((1, 5)).coverage(width=8)) == [(5, 1), (3, 0)]


@pytest.mark.spec("COV.2")
def test_coverage_overlap_counts() -> None:
    assert runs(ir((1, 3), (1, 3), (1, 3)).coverage()) == [(3, 3)]
    assert runs(ir((1, 5), (3, 8)).coverage(width=7)) == [(2, 1), (3, 2), (2, 1)]


@pytest.mark.spec("COV.3")
def test_rle_structure() -> None:
    r = ir((1, 3), (4, 6)).coverage()
    assert runs(r) == [(6, 1)]
    assert len(r) == 6 and r.to_numpy().tolist() == [1] * 6
    assert r == Rle(np.array([6]), np.array([1]))
    assert r != Rle(np.array([3, 3]), np.array([1, 2]))
    assert len(ir().coverage()) == 0
    with pytest.raises(ValueError):
        Rle(np.array([0]), np.array([1]))
    assert len(ir((1, 3), (5, 4), (10, 13)).coverage().to_numpy()) == 13
