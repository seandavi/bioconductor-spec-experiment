"""int32 boundary validation (ADR-0006): start, end and width must all be valid R integers."""

from __future__ import annotations

import pytest

from ranges_e4 import GRanges, IRanges

M = 2**31 - 1
LO = -M  # lowest valid coordinate; -2^31 is R's NA


def test_extreme_valid_ranges_construct() -> None:
    assert IRanges(start=[LO], width=[5]).end.tolist() == [LO + 4]
    assert IRanges(start=[1], width=[M]).end.tolist() == [M]
    assert IRanges(start=[M], width=[1]).end.tolist() == [M]
    assert IRanges(start=[M], width=[0]).end.tolist() == [M - 1]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"start": [M], "width": [2]},  # end past max
        {"start": [LO], "end": [M]},  # inferred width 2^32-1
        {"start": [LO - 1], "width": [1]},  # start below min
        {"end": [LO - 1], "width": [0]},  # end below min
    ],
)
def test_out_of_range_construction_raises(kwargs: dict) -> None:
    with pytest.raises(ValueError):
        IRanges(**kwargs)


def test_zero_width_end_below_min_raises() -> None:
    with pytest.raises(ValueError):
        IRanges(start=[LO], width=[0])  # end would be -2^31 = NA


def test_zero_width_shift_past_max_raises() -> None:
    x = IRanges(start=[M], width=[0])
    with pytest.raises(ValueError):
        x.shift(1)


def test_zero_width_shift_below_min_raises() -> None:
    x = IRanges(start=[LO + 1], width=[0])  # end = LO, valid
    assert x.shift(0).end.tolist() == [LO]
    with pytest.raises(ValueError):
        x.shift(-1)


@pytest.mark.parametrize("op", ["shift", "resize", "flank", "promoters"])
def test_granges_zero_width_boundaries_raise(op: str) -> None:
    g = GRanges(["a"], start=[M], width=[0])
    call = {
        "shift": lambda: g.shift(1),
        "resize": lambda: g.resize(5),
        "flank": lambda: g.flank(5, start=False),
        "promoters": lambda: g.promoters(0, 5),
    }[op]
    with pytest.raises(ValueError):
        call()
