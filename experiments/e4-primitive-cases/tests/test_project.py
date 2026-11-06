"""Projector result-kind inference."""

from __future__ import annotations

from conformance.ops import result_kind
from conformance.project import project
from ranges_e4 import GRanges


def test_empty_coverage_projects_as_rle_list_when_op_says_so() -> None:
    cov = GRanges([], start=[], width=[]).coverage()
    assert cov == {}
    p = project(cov, False, result_kind({"op": "coverage"}))
    assert p == {"kind": "RleList", "names": [], "elements": []}


def test_empty_named_vector_still_projects_as_vector() -> None:
    g = GRanges([], start=[], width=[])
    p = project(g.seqlengths, False, result_kind({"op": "seqlengths"}))
    assert p["kind"] == "vector" and p["values"] == [] and p["names"] == []
