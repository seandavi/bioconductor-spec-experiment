"""Scope audit: every operation required by TARGET.md is a public API."""

import pytest

from ranges_e0 import GRanges, IRanges, Seqinfo

COMMON = [
    "shift", "narrow", "resize", "flank", "promoters",
    "range", "reduce", "gaps", "disjoin",
    "union", "intersect", "setdiff", "punion", "pintersect", "psetdiff",
    "find_overlaps", "count_overlaps", "overlaps_any", "subset_by_overlaps",
    "nearest", "precede", "follow", "distance", "coverage",
    "order", "sort", "rank", "duplicated", "unique",
    "concat", "with_names", "with_mcols",
]
ACCESSORS = ["start", "end", "width", "names", "mcols"]


@pytest.mark.parametrize("name", COMMON + ACCESSORS)
def test_iranges_has(name):
    assert hasattr(IRanges, name)


@pytest.mark.parametrize("name", COMMON + ACCESSORS + [
    "seqnames", "strand", "ranges", "seqinfo", "seqlevels", "seqlengths", "is_circular", "genome", "trim",
])
def test_granges_has(name):
    assert hasattr(GRanges, name)


def test_seqinfo_has():
    assert hasattr(Seqinfo, "merge")
