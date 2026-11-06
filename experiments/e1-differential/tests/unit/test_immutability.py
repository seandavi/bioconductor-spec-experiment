"""ADR-0002: objects are immutable -- no caller-held array or frame may change one (review finding 2)."""

import numpy as np
import polars as pl
import pytest

from ranges_e1 import GRanges, Hits, IRanges, Rle


def test_constructor_detaches_caller_arrays() -> None:
    base = np.array([1, 5, 9], dtype=np.int64)
    x = IRanges(start=base, width=np.array([2, 2, 2]))
    view = base[:]  # read-only views of a writable base used to leak through
    ro = base.view()
    ro.flags.writeable = False
    y = IRanges(start=ro, width=[1, 1, 1])
    base[0] = 100
    assert x.start.tolist() == [1, 5, 9] and y.start.tolist() == [1, 5, 9]
    assert view[0] == 100


def test_accessor_arrays_cannot_be_reenabled() -> None:
    x = IRanges(start=[1, 2], width=[1, 1])
    for arr in (x.start, x.width):
        with pytest.raises(ValueError):
            arr[0] = 7
        with pytest.raises(ValueError):
            arr.flags.writeable = True
    assert x.end.tolist() == [1, 2]
    g = GRanges("c", start=[1, 2], width=[1, 1])
    with pytest.raises(ValueError):
        g.start.flags.writeable = True


def test_mcols_frames_are_isolated_both_ways() -> None:
    df = pl.DataFrame({"m": [1]})
    x = IRanges(start=[1], width=[1], mcols=df)
    shifted = x.shift(1)
    df.vstack(pl.DataFrame({"m": [2]}), in_place=True)  # caller mutates its own frame
    assert x.mcols is not None and x.mcols.height == 1 and shifted.mcols is not None and shifted.mcols.height == 1
    got = x.mcols
    assert got is not None
    got.vstack(pl.DataFrame({"m": [3]}), in_place=True)  # caller mutates what the accessor returned
    assert x.mcols is not None and x.mcols.height == 1 and shifted.mcols.height == 1
    g = GRanges("c", start=[1], width=[1], mcols=df.head(1))
    g_df = g.mcols
    assert g_df is not None
    g_df.vstack(pl.DataFrame({"m": [9]}), in_place=True)
    assert g.mcols is not None and g.mcols.height == 1


def test_hits_and_rle_do_not_freeze_or_alias_caller_arrays() -> None:
    q, s = np.array([0, 1]), np.array([1, 0])
    h = Hits(q, s, 2, 2)
    q[0] = 5
    assert h.query_hits.tolist() == [0, 1] and q.flags.writeable
    v, n = np.array([1, 2]), np.array([3, 4])
    r = Rle(v, n)
    v[0] = 9
    assert r.values.tolist() == [1, 2] and v.flags.writeable
