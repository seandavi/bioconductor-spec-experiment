"""Python behaviour where R is skipped as defective (ADR-0010) or Python deliberately raises."""

from __future__ import annotations

import numpy as np
import pytest

from ranges_e4 import GRanges, IRanges, Seqinfo

M = 2**31 - 1


def test_disjoin_with_end_at_int_max() -> None:
    d = IRanges(start=[1, M - 5], width=[M, 6]).disjoin()
    assert (d.start.tolist(), d.end.tolist()) == ([1], [M])


def test_find_overlaps_with_maxgap_int_max_hits() -> None:
    h = IRanges(start=[1], width=[1]).find_overlaps(IRanges(start=[M - 5], width=[1]), maxgap=M)
    assert (h.query.tolist(), h.subject.tolist()) == ([0], [0])


def test_gaps_and_setdiff_at_int_max() -> None:
    g = IRanges(start=[5], width=[3]).gaps(start=1, end=M)
    assert (g.start.tolist(), g.end.tolist()) == ([1, 8], [4, M])
    s = IRanges(start=[1], width=[M]).setdiff(IRanges(start=[5], width=[3]))
    assert (s.start.tolist(), s.end.tolist()) == ([1, 8], [4, M])


def test_empty_query_precede_follow_return_select_shaped_vectors() -> None:
    e = GRanges([], start=[], width=[], seqinfo=Seqinfo(["c"], [30]))
    g = GRanges(["c"], start=[5], width=[3], seqinfo=Seqinfo(["c"], [30]))
    for fn in (e.precede, e.follow):
        r = fn(g)
        assert isinstance(r, np.ma.MaskedArray) and len(r) == 0


@pytest.mark.parametrize("op", ["resize", "flank", "promoters"])
def test_results_beyond_int32_raise(op: str) -> None:
    x = IRanges(start=[M - 5], width=[3])
    with pytest.raises(ValueError):
        {"resize": lambda: x.resize(100), "flank": lambda: x.flank(100, start=False), "promoters": lambda: x.promoters(0, 100)}[op]()


def test_distance_beyond_int32_is_masked() -> None:
    d = IRanges(start=[-M + 5], width=[1]).distance(IRanges(start=[M - 5], width=[1]))
    assert d.mask.tolist() == [True]


def test_circular_coverage_wraps() -> None:
    g = GRanges(["c"], start=[25], width=[10], seqinfo=Seqinfo(["c"], [30], [True]))
    cov = g.coverage()["c"]
    assert cov.values.tolist() == [1, 0, 1] and cov.lengths.tolist() == [4, 20, 6]
