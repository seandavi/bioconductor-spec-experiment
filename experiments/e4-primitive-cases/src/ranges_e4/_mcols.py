"""Per-range metadata columns, stored as a polars DataFrame.

A frame with no columns has height 0 whatever the object length, so
"no mcols" is represented by an empty frame and never row-checked.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import polars as pl


def empty() -> pl.DataFrame:
    return pl.DataFrame()


def coerce(mcols: pl.DataFrame | Mapping[str, Sequence[Any]] | None, n: int) -> pl.DataFrame:
    """Accept a frame or a dict of columns; check height against ``n``."""
    if mcols is None:
        return empty()
    df = mcols if isinstance(mcols, pl.DataFrame) else pl.DataFrame(dict(mcols))
    if df.width and df.height != n:
        raise ValueError(f"mcols has {df.height} rows, expected {n}")
    return df


def take(df: pl.DataFrame, idx: np.ndarray) -> pl.DataFrame:
    if not df.width:
        return df
    return df[idx.astype(np.int64)]


def concat(frames: Sequence[pl.DataFrame], lengths: Sequence[int]) -> pl.DataFrame:
    """Row-bind; columns missing from a frame become null (R fills NA)."""
    cols: list[str] = []
    for f in frames:
        cols.extend(c for c in f.columns if c not in cols)
    if not cols:
        return empty()
    parts = [
        pl.DataFrame({c: f[c] if c in f.columns else pl.Series(c, [None] * n) for c in cols})
        for f, n in zip(frames, lengths, strict=True)
    ]
    return pl.concat(parts, how="vertical_relaxed")
