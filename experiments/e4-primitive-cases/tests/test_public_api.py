"""Every operation required by TARGET.md is a public attribute of the classes."""

from __future__ import annotations

import pytest

from ranges_e4 import GRanges, Hits, IRanges, Rle, Seqinfo

BOTH = [
    "start", "end", "width", "names", "mcols", "concat", "shift", "narrow", "resize", "flank",
    "promoters", "range", "reduce", "gaps", "disjoin", "union", "intersect", "setdiff", "punion",
    "pintersect", "psetdiff", "find_overlaps", "count_overlaps", "overlaps_any",
    "subset_by_overlaps", "nearest", "precede", "follow", "distance", "coverage", "order", "sort",
    "rank", "duplicated", "unique", "with_mcols",
]
GR_ONLY = ["seqnames", "strand", "ranges", "seqinfo", "seqlevels", "seqlengths", "is_circular", "genome", "trim"]


@pytest.mark.parametrize("cls", [IRanges, GRanges])
@pytest.mark.parametrize("name", BOTH)
def test_shared_operation_is_public(cls: type, name: str) -> None:
    assert hasattr(cls, name)


@pytest.mark.parametrize("name", GR_ONLY)
def test_granges_operation_is_public(name: str) -> None:
    assert hasattr(GRanges, name)


def test_supporting_types_exported() -> None:
    assert all(t is not None for t in (Hits, Rle, Seqinfo))
    for cls in (IRanges, GRanges):
        assert hasattr(cls, "__len__") and hasattr(cls, "__getitem__")
