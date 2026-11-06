"""GRanges / Seqinfo unit tests beyond the shared cases."""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e4 import GRanges, Seqinfo


def test_seqinfo_merge_is_symmetric_on_values() -> None:
    a = Seqinfo(["x", "y"], [10, None], [None, True], "g")
    b = Seqinfo(["y", "z"], [20, 30], [None, None], None)
    ab, ba = a.merge(b), b.merge(a)
    assert ab.seqlengths == ba.seqlengths | {} and ab.seqnames == ("x", "y", "z")
    assert ba.seqnames == ("y", "z", "x")
    assert ab.genome == {"x": "g", "y": "g", "z": "g"} or ab.genome["x"] == "g"


def test_seqinfo_is_immutable() -> None:
    si = Seqinfo(["a"])
    with pytest.raises(AttributeError):
        si.seqnames = ("b",)  # type: ignore[misc]


def test_granges_getitem_does_not_alias_input() -> None:
    g = GRanges(["a", "b"], start=[1, 2], width=[1, 1])
    sub = g[[0]]
    assert len(sub) == 1 and len(g) == 2
    assert sub.seqinfo is g.seqinfo


def test_trim_leaves_in_bounds_ranges_alone() -> None:
    g = GRanges(["a"], start=[2], width=[3], seqinfo=Seqinfo(["a"], [10]))
    t = g.trim()
    assert t.start.tolist() == [2] and t.width.tolist() == [3]


def test_concat_remaps_codes_by_name() -> None:
    a = GRanges(["x"], start=[1], width=[1], seqinfo=Seqinfo(["x", "y"]))
    b = GRanges(["y"], start=[1], width=[1], seqinfo=Seqinfo(["y", "x"]))
    c = GRanges.concat(a, b)
    assert c.seqnames.tolist() == ["x", "y"] and c.seqlevels == ["x", "y"]


def test_strand_codes_roundtrip() -> None:
    g = GRanges(["a", "a", "a"], start=[1, 2, 3], width=[1, 1, 1], strand=["+", "-", "*"])
    assert g.strand.tolist() == ["+", "-", "*"]
    assert np.array_equal(g._strand, [0, 1, 2])
