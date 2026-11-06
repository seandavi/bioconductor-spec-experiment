"""Row-aligned metadata columns (polars DataFrame) helpers.

polars cannot hold zero columns with non-zero height, so "no columns" is represented by
`None` in the owning object; `n` (the object length) is always the row count.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from ._util import IntArray


def normalise(df: pl.DataFrame | None, n: int) -> pl.DataFrame | None:
    """Validate/recycle a user-supplied frame to n rows; zero columns -> None.

    R recycles a shorter frame cyclically; more rows than the object (or 0 rows for n > 0) is an error.
    """
    if df is None or df.width == 0:
        return None
    h = df.height
    if h == n:
        return df.clone()  # a private wrapper: in-place edits (replace_column) must not leak
    if h == 0 or h > n:
        raise ValueError(f"trying to set a metadata column of length {h} on an object of length {n}")
    return df[np.resize(np.arange(h), n)]


def take(df: pl.DataFrame | None, idx: IntArray) -> pl.DataFrame | None:
    if df is None:
        return None
    return df[idx]


def concat(parts: list[tuple[int, pl.DataFrame | None]]) -> pl.DataFrame | None:
    """Row-bind frames; columns are unioned in order of first appearance, gaps filled with NA.

    Column types are widened to a common supertype (integer + double -> double); incompatible
    types raise rather than truncate.
    """
    order: list[str] = []
    for _, df in parts:
        if df is not None:
            order.extend(c for c in df.columns if c not in order)
    if not order:
        return None
    frames = []
    for n, df in parts:
        cols = [
            df[name] if df is not None and name in df.columns else pl.Series(name, [None] * n, dtype=pl.Null)
            for name in order
        ]
        frames.append(pl.DataFrame(cols))
    try:
        return pl.concat(frames, how="vertical_relaxed")
    except (pl.exceptions.PolarsError, TypeError) as e:
        raise ValueError(f"cannot combine metadata columns with incompatible types: {e}") from None


def empty_frame() -> pl.DataFrame:
    return pl.DataFrame()
