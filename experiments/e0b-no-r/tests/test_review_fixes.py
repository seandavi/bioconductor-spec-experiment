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


