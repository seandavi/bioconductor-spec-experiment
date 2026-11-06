"""Regressions from the independent review (findings 1-4)."""

import numpy as np
import polars as pl
import pytest

from ranges_e0b import GRanges, IRanges, Seqinfo


# --- finding 1: ignore_strand on strand-aware intra ops ---------------------------


def minus_plus() -> GRanges:
    return GRanges("c", IRanges(start=[100, 100], end=[110, 110]), ["-", "+"])


def coords(x: GRanges) -> list[tuple[int, int]]:
    return list(zip(x.start.tolist(), x.end.tolist()))


def test_resize_ignore_strand_matches_iranges() -> None:
    x = minus_plus()
    assert coords(x.resize(3)) == [(108, 110), (100, 102)]
    assert coords(x.resize(3, ignore_strand=True)) == coords(GRanges("c", x.ranges).resize(3))
    assert coords(x.resize(3, fix="end", ignore_strand=True)) == [(108, 110)] * 2


def test_flank_ignore_strand_matches_iranges() -> None:
    x = minus_plus()
    assert coords(x.flank(5)) == [(111, 115), (95, 99)]
    assert coords(x.flank(5, ignore_strand=True)) == [(95, 99)] * 2
    assert coords(x.flank(5, start=False, both=True, ignore_strand=True)) == [(106, 115)] * 2


def test_promoters_ignore_strand_matches_iranges() -> None:
    x = minus_plus()
    assert coords(x.promoters(10, 5)) == [(106, 120), (90, 104)]
    assert coords(x.promoters(10, 5, ignore_strand=True)) == [(90, 104)] * 2


# --- finding 2: mcols ownership ---------------------------------------------------


@pytest.mark.parametrize("make", [
    lambda m: IRanges(start=[1], width=2, mcols=m),
    lambda m: GRanges("c", IRanges(start=[1], width=2), mcols=m),
])
def test_mcols_not_aliased_to_caller_or_accessor(make) -> None:  # type: ignore[no-untyped-def]
    m = pl.DataFrame({"score": [1]})
    x = make(m)
    y = x.shift(5)
    m.replace_column(0, pl.Series("score", [99]))  # caller mutates their frame
    assert x.mcols["score"].to_list() == [1] and y.mcols["score"].to_list() == [1]
    x.mcols.insert_column(1, pl.Series("extra", [0]))  # mutating the accessor result
    x.mcols.replace_column(0, pl.Series("score", [7]))
    assert x.mcols.columns == ["score"]
    assert x.mcols["score"].to_list() == [1] and y.mcols["score"].to_list() == [1]


def test_with_mcols_copies() -> None:
    m = pl.DataFrame({"v": [1, 2]})
    x = IRanges(start=[1, 2], width=1).with_mcols(m)
    m.replace_column(0, pl.Series("v", [0, 0]))
    assert x.mcols["v"].to_list() == [1, 2]


