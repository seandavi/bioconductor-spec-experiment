"""IRanges: integer ranges, 1-based and closed (end = start + width - 1)."""

from collections.abc import Mapping, Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from ._sew import (
    IntArray,
    concat_mcols,
    concat_names,
    normalize_mcols,
    resolve_index,
    solve_sew,
    take_mcols,
)


class IRanges:
    """Immutable vector of integer ranges with optional names and mcols.

    Built from any two of `start`, `end`, `width`. Positional indices in all
    results are 0-based; coordinates are 1-based like R.
    """

    __slots__ = ("_start", "_width", "_names", "_mcols")

    def __init__(
        self,
        start: Any = None,
        end: Any = None,
        width: Any = None,
        names: Sequence[str] | None = None,
        mcols: pl.DataFrame | Mapping[str, Any] | None = None,
    ) -> None:
        s, w = solve_sew(start, end, width)
        self._init(s.copy(), w.copy(), names, mcols)

    def _init(
        self,
        start: IntArray,
        width: IntArray,
        names: Sequence[str] | None,
        mcols: pl.DataFrame | Mapping[str, Any] | None,
    ) -> None:
        n = len(start)
        start.setflags(write=False)
        width.setflags(write=False)
        self._start = start
        self._width = width
        if names is not None:
            names = list(names)
            if len(names) != n:
                raise ValueError(f"names has length {len(names)}, expected {n}")
        self._names = names
        self._mcols = normalize_mcols(mcols, n)

    @classmethod
    def _new(
        cls,
        start: IntArray,
        width: IntArray,
        names: list[str] | None = None,
        mcols: pl.DataFrame | None = None,
    ) -> Self:
        """Internal constructor from validated arrays (width >= 0)."""
        obj = object.__new__(cls)
        obj._init(np.ascontiguousarray(start, dtype=np.int64), np.ascontiguousarray(width, dtype=np.int64), names, mcols)
        return obj

    # --- accessors -----------------------------------------------------------

    @property
    def start(self) -> IntArray:
        return self._start

    @property
    def width(self) -> IntArray:
        return self._width

    @property
    def end(self) -> IntArray:
        return self._start + self._width - 1

    @property
    def names(self) -> list[str] | None:
        return None if self._names is None else list(self._names)

    @property
    def mcols(self) -> pl.DataFrame:
        return self._mcols

    def __len__(self) -> int:
        return len(self._start)

    def __repr__(self) -> str:
        n = len(self)
        head = f"IRanges of length {n}"
        if self._mcols.width:
            head += f" with {self._mcols.width} mcols ({', '.join(self._mcols.columns)})"
        rows = [f"{head}"]
        shown = range(n) if n <= 10 else list(range(5)) + [-1] + list(range(n - 5, n))
        for i in shown:
            if i == -1:
                rows.append("  ...")
                continue
            nm = f" {self._names[i]}" if self._names is not None else ""
            rows.append(f"  [{i}] {self._start[i]}..{self._start[i] + self._width[i] - 1} (w={self._width[i]}){nm}")
        return "\n".join(rows)

    def with_names(self, names: Sequence[str] | None) -> Self:
        return self._new(self._start, self._width, None if names is None else list(names), self._mcols)

    def with_mcols(self, mcols: pl.DataFrame | Mapping[str, Any] | None) -> Self:
        return self._new(self._start, self._width, self._names, normalize_mcols(mcols, len(self)))

    def _replace(self, start: IntArray, width: IntArray) -> Self:
        """Same names/mcols, new coordinates (length must be unchanged)."""
        return self._new(start, width, self._names, self._mcols)

    def equals(self, other: "IRanges") -> bool:
        """Same starts, widths, names and mcols (there is no element-wise `==`)."""
        return (
            len(self) == len(other)
            and bool(np.array_equal(self._start, other._start))
            and bool(np.array_equal(self._width, other._width))
            and self._names == other._names
            and self._mcols.equals(other._mcols)
        )

    # --- subsetting and combining -------------------------------------------

    def __getitem__(self, i: Any) -> Self:
        idx = resolve_index(i, len(self), self._names)
        return self._take(idx)

    def _take(self, idx: IntArray) -> Self:
        names = None if self._names is None else [self._names[j] for j in idx]
        return self._new(self._start[idx], self._width[idx], names, take_mcols(self._mcols, idx))

    @classmethod
    def combine(cls, items: Sequence["IRanges"]) -> Self:
        """Concatenate; mcols must have identical column names."""
        if not items:
            return cls._new(np.empty(0, np.int64), np.empty(0, np.int64))
        return cls._new(
            np.concatenate([x._start for x in items]),
            np.concatenate([x._width for x in items]),
            concat_names([x._names for x in items], [len(x) for x in items]),
            concat_mcols([x._mcols for x in items]),
        )

    def concat(self, *others: "IRanges") -> Self:
        return self.combine([self, *others])
