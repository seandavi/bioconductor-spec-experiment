"""Cross-cutting conventions (spec/00-conventions.md); per-op cases live with their area."""

from __future__ import annotations

import polars as pl
import pytest

import ranges_e3
from helpers import ir, se


@pytest.mark.spec("CONV.4")
def test_missing_ints_are_masked_arrays() -> None:
    import numpy as np

    m = ir((1, 3)).find_overlaps(ir((9, 9)), select="first")
    assert isinstance(m, np.ma.MaskedArray) and m.mask.tolist() == [True]


@pytest.mark.spec("CONV.5")
def test_inputs_never_mutated() -> None:
    x = ir((5, 6), (1, 2))
    before = se(x)
    x.shift(3)
    x.reduce()
    x.sort()
    assert se(x) == before


@pytest.mark.spec("CONV.6")
def test_per_range_argument_length_rule() -> None:
    x = ir((1, 3), (5, 6))
    assert se(x.shift([1, 2])) == [(2, 4), (7, 8)]
    assert se(x.shift([1])) == [(2, 4), (6, 7)]
    with pytest.raises(ValueError):
        x.shift([1, 2, 3])


@pytest.mark.spec("CONV.8")
def test_names_mcols_propagation_table() -> None:
    x = ir((1, 3), (5, 6), names=["a", "b"]).with_mcols(pl.DataFrame({"k": [1, 2]}))
    for y in (x.shift(1), x.narrow(start=1), x.resize(1), x.flank(1), x.promoters(1, 1), x.sort(), x[[1, 0]]):
        assert y.names is not None and y.mcols.columns == ["k"]
    for y in (x.reduce(), x.disjoin(), x.gaps(start=1, end=10), x.range(), x.union(x)):
        assert y.names is None and y.mcols.width == 0
    p = x.punion(x)
    assert p.names == ["a", "b"] and p.mcols.width == 0


@pytest.mark.spec("CONV.9")
def test_public_names() -> None:
    for name in ("find_overlaps", "count_overlaps", "overlaps_any", "subset_by_overlaps"):
        assert hasattr(ranges_e3.IRanges, name)
    assert callable(ranges_e3.concat)
