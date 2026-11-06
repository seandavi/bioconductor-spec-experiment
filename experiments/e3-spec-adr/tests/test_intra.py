"""Intra-range operations (spec/20-intra.md)."""

from __future__ import annotations

import polars as pl
import pytest

from helpers import ir, se


@pytest.mark.spec("INTRA.1")
def test_shift() -> None:
    assert se(ir((1, 3), (5, 4)).shift(2)) == [(3, 5), (7, 6)]
    assert se(ir((1, 3), (5, 4), (10, 13), (20, 20)).shift([1, -1, 0, 10])) == [
        (2, 4), (4, 3), (10, 13), (30, 30),
    ]


@pytest.mark.spec("INTRA.2")
def test_narrow_relative_positions() -> None:
    x = ir((10, 13))
    assert se(x.narrow(start=2)) == [(11, 13)]
    assert se(x.narrow(start=-2)) == [(12, 13)]
    assert se(x.narrow(end=-2)) == [(10, 12)]
    assert se(x.narrow(end=2)) == [(10, 11)]
    assert se(x.narrow(start=2, end=3)) == [(11, 12)]
    assert se(x.narrow(start=2, width=2)) == [(11, 12)]
    assert se(x.narrow(end=3, width=2)) == [(11, 12)]
    assert se(x.narrow(start=-1)) == [(13, 13)]
    assert se(x.narrow(start=-3, end=-2)) == [(11, 12)]


@pytest.mark.spec("INTRA.2")
def test_narrow_zero_width_results() -> None:
    assert se(ir((1, 3)).narrow(start=3, width=0)) == [(3, 2)]
    assert se(ir((10, 13)).narrow(end=0)) == [(10, 9)]
    assert se(ir((10, 13)).narrow(start=5)) == [(14, 13)]
    assert se(ir((10, 13)).narrow(end=-5)) == [(10, 9)]
    assert se(ir((5, 4)).narrow(start=1)) == [(5, 4)]


@pytest.mark.spec("INTRA.2")
def test_narrow_without_arguments_is_identity() -> None:
    assert se(ir((10, 13), (5, 4)).narrow()) == [(10, 13), (5, 4)]


@pytest.mark.spec("INTRA.2")
@pytest.mark.parametrize(
    "kwargs",
    [{"width": 3}, {"start": 1, "end": 3, "width": 3}, {"start": -5}, {"start": 6},
     {"start": 4, "end": 2}, {"end": 5}, {"end": -6}, {"start": 1, "width": 5}, {"end": 4, "width": 5}],
)
def test_narrow_refusals(kwargs: dict[str, int]) -> None:
    with pytest.raises(ValueError):
        ir((10, 13)).narrow(**kwargs)
    # zero-width input only admits start=1
    with pytest.raises(ValueError):
        ir((5, 4)).narrow(start=2)


@pytest.mark.spec("INTRA.3")
def test_resize_fix_modes() -> None:
    x = ir((1, 3), (5, 4), (10, 13), (20, 20))
    assert se(x.resize(2)) == [(1, 2), (5, 6), (10, 11), (20, 21)]
    assert se(x.resize(2, "end")) == [(2, 3), (3, 4), (12, 13), (19, 20)]
    assert se(x.resize(2, "center")) == [(1, 2), (4, 5), (11, 12), (19, 20)]
    assert se(x.resize(0, "end")) == [(4, 3), (5, 4), (14, 13), (21, 20)]
    assert se(x.resize(5, "center")) == [(0, 4), (2, 6), (9, 13), (18, 22)]
    assert se(x.resize([1, 2, 3, 4])) == [(1, 1), (5, 6), (10, 12), (20, 23)]


@pytest.mark.spec("INTRA.3")
@pytest.mark.parametrize(
    "start,end,new,expected",
    [(1, 5, 2, (2, 3)), (1, 4, 3, (1, 3)), (1, 4, 2, (2, 3)), (1, 3, 2, (1, 2)),
     (2, 4, 0, (3, 2)), (1, 6, 3, (2, 4)), (1, 6, 1, (3, 3)), (1, 5, 4, (1, 4))],
)
def test_resize_center_rounding(start: int, end: int, new: int, expected: tuple[int, int]) -> None:
    assert se(ir((start, end)).resize(new, "center")) == [expected]


@pytest.mark.spec("INTRA.3")
def test_resize_refusals() -> None:
    with pytest.raises(ValueError):
        ir((1, 3)).resize(-1)
    with pytest.raises(ValueError):
        ir((1, 3)).resize(1, "middle")


@pytest.mark.spec("INTRA.4")
def test_flank_one_sided() -> None:
    x = ir((1, 3), (5, 4), (10, 13), (20, 20))
    assert se(x.flank(2)) == [(-1, 0), (3, 4), (8, 9), (18, 19)]
    assert se(x.flank(-2)) == [(1, 2), (5, 6), (10, 11), (20, 21)]
    assert se(x.flank(0)) == [(1, 0), (5, 4), (10, 9), (20, 19)]
    assert se(x.flank(2, start=False)) == [(4, 5), (5, 6), (14, 15), (21, 22)]
    assert se(x.flank(-2, start=False)) == [(2, 3), (3, 4), (12, 13), (19, 20)]
    assert se(x.flank(0, start=False)) == [(4, 3), (5, 4), (14, 13), (21, 20)]


@pytest.mark.spec("INTRA.4")
def test_flank_both_and_per_range_side() -> None:
    x = ir((1, 3), (5, 4), (10, 13), (20, 20))
    assert se(x.flank(2, both=True)) == [(-1, 2), (3, 6), (8, 11), (18, 21)]
    assert se(x.flank(2, start=False, both=True)) == [(2, 5), (3, 6), (12, 15), (19, 22)]
    assert se(x.flank(-2, both=True)) == se(x.flank(2, both=True))
    assert se(x.flank(1, both=True)) == [(0, 1), (4, 5), (9, 10), (19, 20)]
    assert se(ir((10, 14), (10, 14)).flank([2, 3], start=[True, False])) == [(8, 9), (15, 17)]


@pytest.mark.spec("INTRA.5")
def test_promoters() -> None:
    x = ir((1, 3), (5, 4))
    assert se(x.promoters(3, 1)) == [(-2, 1), (2, 5)]
    assert se(x.promoters()) == [(-1999, 200), (-1995, 204)]
    assert se(x.promoters(0, 0)) == [(1, 0), (5, 4)]
    with pytest.raises(ValueError):
        x.promoters(-1, 0)


@pytest.mark.spec("CONV.8")
def test_intra_keeps_names_and_mcols() -> None:
    x = ir((1, 3), (5, 6), names=["a", "b"]).with_mcols(pl.DataFrame({"k": [1, 2]}))
    for y in (x.shift(1), x.narrow(start=1), x.resize(1), x.flank(1), x.promoters(1, 1)):
        assert y.names == ["a", "b"] and y.mcols["k"].to_list() == [1, 2]
