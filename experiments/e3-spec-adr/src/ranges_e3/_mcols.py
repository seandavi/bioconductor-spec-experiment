"""mcols helpers (ADR-0003): ``None`` internally means "no columns"."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import polars as pl

from ranges_e3._types import IntArray


def validate(mcols: pl.DataFrame | None, n: int) -> pl.DataFrame | None:
    """Normalise a user-supplied table: no columns -> ``None``; wrong height refused (CORE.8)."""
    if mcols is None or mcols.width == 0:
        return None
    if mcols.height != n:
        raise ValueError(f"mcols has {mcols.height} rows, expected {n}")
    return mcols


def public(mcols: pl.DataFrame | None) -> pl.DataFrame:
    """The table shown to users: empty frame when there are no columns."""
    return pl.DataFrame() if mcols is None else mcols


def take(mcols: pl.DataFrame | None, index: IntArray) -> pl.DataFrame | None:
    """Rows ``index`` of ``mcols`` (CORE.10)."""
    if mcols is None:
        return None
    return mcols[index.tolist()] if len(index) else mcols.head(0)


def concat(parts: Sequence[tuple[pl.DataFrame | None, int]]) -> pl.DataFrame | None:
    """Stack tables over the union of columns, null-filling gaps (CORE.12)."""
    if all(m is None for m, _ in parts):
        return None
    columns: list[str] = []
    for m, _ in parts:
        if m is not None:
            columns.extend(c for c in m.columns if c not in columns)
    frames = []
    for m, n in parts:
        if m is None:
            frames.append(pl.DataFrame({c: [None] * n for c in columns}))
        else:
            frames.append(m)
    return pl.concat(frames, how="diagonal_relaxed").select(columns)


def with_column(mcols: pl.DataFrame | None, name: str, values: pl.Series) -> pl.DataFrame:
    """Add or replace one column."""
    base = pl.DataFrame() if mcols is None else mcols
    return base.with_columns(values.alias(name))


def revmap_series(groups: Sequence[IntArray]) -> pl.Series:
    """List(Int64) column from per-range index arrays (INTER.5)."""
    return pl.Series("revmap", [g.tolist() for g in groups], dtype=pl.List(pl.Int64))


def empty_int() -> IntArray:
    return np.empty(0, dtype=np.int64)
