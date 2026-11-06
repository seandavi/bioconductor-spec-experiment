"""Comparator regression tests: exact means exact; set/multiset do not hide other differences."""

from __future__ import annotations

from typing import Any

from conformance.compare import compare


def vec(values: list[Any]) -> dict[str, Any]:
    return {"kind": "vector", "values": values, "names": None}


def ir(start: list[int], width: list[int], **extra: Any) -> dict[str, Any]:
    end = [s + w - 1 for s, w in zip(start, width, strict=True)]
    return {"kind": "IRanges", "start": start, "end": end, "width": width, "names": None, "mcols": {}} | extra


def test_exact_rejects_off_by_one_integer() -> None:
    assert compare({}, vec([2147483647]), vec([2147483646]))
    assert compare({}, vec([1, 2]), vec([1, 3]))


def test_exact_rejects_unequal_floats() -> None:
    assert compare({}, vec([0.5]), vec([0.5000001]))
    assert not compare({}, vec([0.5]), vec([0.5]))


def test_exact_int_equals_same_valued_float() -> None:
    assert not compare({}, vec([3]), vec([3.0]))


def test_exact_distinguishes_bool_from_int_and_na() -> None:
    assert compare({}, vec([True]), vec([1]))
    assert compare({}, vec([None]), vec([0]))


def test_set_mode_ignores_member_order_only() -> None:
    spec = {"compare": {"fields": {"mcols.revmap": "set"}}}
    g = ir([1, 5], [1, 1], mcols={"revmap": [[1, 2], [3]]})
    same = ir([1, 5], [1, 1], mcols={"revmap": [[2, 1], [3]]})
    swapped_rows = ir([1, 5], [1, 1], mcols={"revmap": [[3], [1, 2]]})
    assert not compare(spec, g, same)
    assert compare(spec, g, swapped_rows)


def test_multiset_ignores_row_order() -> None:
    spec = {"compare": {"rows": "multiset"}}
    assert not compare(spec, ir([1, 5], [1, 2]), ir([5, 1], [2, 1]))
    assert compare(spec, ir([1, 5], [1, 2]), ir([5, 1], [2, 2]))


def test_multiset_rejects_truncated_rows() -> None:
    spec = {"compare": {"rows": "multiset"}}
    bad = {"kind": "IRanges", "start": [1, 99], "end": [1], "width": [1], "names": None, "mcols": {}}
    good = {"kind": "IRanges", "start": [1], "end": [1], "width": [1], "names": None, "mcols": {}}
    assert compare(spec, good, bad)
    assert compare(spec, bad, good)
