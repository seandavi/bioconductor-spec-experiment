"""Unit tests for intra-range ops beyond the shared cases."""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e4 import IRanges


def test_resize_center_matches_floor_rule() -> None:
    x = IRanges(start=[5, 5], width=[4, 3])
    r = x.resize(2, fix="center")
    assert r.start.tolist() == [6, 5]


def test_flank_both_is_union_of_inner_and_outer() -> None:
    x = IRanges(start=[10], width=[5])
    both = x.flank(3, both=True)
    outer = x.flank(3)
    inner = x.flank(-3)
    assert both.start[0] == outer.start[0]
    assert both.end[0] == inner.end[0]


def test_shift_rejects_wrong_length_vector() -> None:
    with pytest.raises(ValueError):
        IRanges(start=[1, 2, 3], width=[1, 1, 1]).shift([1, 2])


def test_intra_ops_do_not_mutate_input() -> None:
    x = IRanges(start=[1, 2], width=[3, 4])
    before = (x.start.copy(), x.width.copy())
    x.shift(5).resize(1).narrow(start=1, width=1).flank(2).promoters(1, 1)
    assert (x.start == before[0]).all() and (x.width == before[1]).all()


def test_non_integer_float_rejected() -> None:
    with pytest.raises(ValueError):
        IRanges(start=np.array([1.5]), width=[1])
