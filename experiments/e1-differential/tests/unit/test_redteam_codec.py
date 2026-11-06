"""Regression tests for typed empty metadata in differential/codec.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "differential"))
import adapter  # noqa: E402
import codec  # noqa: E402


@pytest.mark.parametrize("typ,dtype", [("int", pl.Int64), ("num", pl.Float64), ("lgl", pl.Boolean), ("chr", pl.String)])
def test_codec_empty_metadata_keeps_declared_type(typ: str, dtype: pl.DataType) -> None:
    frame = codec.decode({"DataFrame": {"score": {typ: []}}})
    assert frame.schema["score"] == dtype
    assert adapter.encode(frame)["columns"]["score"]["type"] == typ


@pytest.mark.parametrize("dtype", [np.int64, np.float64, np.bool_])
def test_codec_empty_masked_metadata_keeps_dtype(dtype: type[np.generic]) -> None:
    empty = np.ma.masked_array(np.array([], dtype=dtype), mask=[])
    frame = codec._series("score", empty).to_frame()
    assert frame.schema["score"] != pl.Null
