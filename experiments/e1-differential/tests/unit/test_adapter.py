"""Unit tests for the differential adapter's index translation hooks."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "differential"))
import adapter  # noqa: E402


def test_r_subscript_positive_is_zero_based() -> None:
    assert adapter.r_subscript(np.array([3, 1]), 3).tolist() == [2, 0]


def test_r_subscript_zero_dropped() -> None:
    assert adapter.r_subscript(np.array([0, 2]), 3).tolist() == [1]


def test_r_subscript_negative_excludes() -> None:
    assert adapter.r_subscript(np.array([-1, -3]), 3).tolist() == [False, True, False]


def test_r_subscript_logical_and_names_pass_through() -> None:
    assert adapter.r_subscript(np.array([True, False]), 2).tolist() == [True, False]
    assert adapter.r_subscript(["a"], 2) == ["a"]


def test_to_one_based_keeps_mask() -> None:
    m = np.ma.masked_array([0, 5], mask=[False, True])
    out = adapter.to_one_based(m)
    assert out[0] == 1 and out.mask.tolist() == [False, True]


def test_snake() -> None:
    assert adapter.snake("findOverlaps") == "find_overlaps"
    assert adapter.snake("min.gapwidth") == "min_gapwidth"
    assert adapter.snake("ignore.strand") == "ignore_strand"
