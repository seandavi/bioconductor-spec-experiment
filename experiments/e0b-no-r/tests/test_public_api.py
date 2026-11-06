import pytest

import ranges_e0b
from ranges_e0b import GRanges, IRanges

SHARED = [
    "start", "end", "width", "names", "mcols", "shift", "narrow", "resize", "flank", "promoters",
    "range", "reduce", "gaps", "disjoin", "union", "intersect", "setdiff", "punion", "pintersect",
    "psetdiff", "find_overlaps", "count_overlaps", "overlaps_any", "subset_by_overlaps", "nearest",
    "precede", "follow", "distance", "coverage", "order", "sort", "rank", "duplicated", "unique",
    "concat", "combine", "equals", "with_names", "with_mcols", "__len__", "__getitem__",
]
GRANGES_ONLY = ["seqnames", "strand", "ranges", "seqinfo", "seqlevels", "seqlengths", "is_circular",
                "genome", "trim", "with_seqinfo"]


@pytest.mark.parametrize("cls", [IRanges, GRanges])
def test_required_operations_exist(cls: type) -> None:
    missing = [n for n in SHARED if not hasattr(cls, n)]
    assert not missing


def test_granges_specific_operations_exist() -> None:
    assert not [n for n in GRANGES_ONLY if not hasattr(GRanges, n)]


def test_package_exports() -> None:
    for name in ranges_e0b.__all__:
        assert hasattr(ranges_e0b, name)
    assert {"IRanges", "GRanges", "Seqinfo", "Hits", "Rle", "NO_HIT"} <= set(ranges_e0b.__all__)


STRAND_AWARE = ["resize", "flank", "promoters", "range", "reduce", "gaps", "disjoin", "union", "intersect",
                "setdiff", "punion", "pintersect", "psetdiff", "find_overlaps", "count_overlaps",
                "overlaps_any", "subset_by_overlaps", "nearest", "precede", "follow", "distance"]


@pytest.mark.parametrize("name", STRAND_AWARE)
def test_granges_strand_aware_ops_accept_ignore_strand(name: str) -> None:
    import inspect

    assert "ignore_strand" in inspect.signature(getattr(GRanges, name)).parameters
