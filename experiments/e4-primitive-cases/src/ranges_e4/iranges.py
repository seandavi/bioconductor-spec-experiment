"""IRanges: 1-based, closed integer ranges (``end = start + width - 1``)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _mcols
from ._arrays import as_int, recycle


def as_names(names: Sequence[str | None] | np.ndarray | None, n: int) -> np.ndarray | None:
    if names is None:
        return None
    arr = np.empty(len(names), dtype=object)
    arr[:] = list(names)
    if len(arr) != n:
        raise ValueError(f"names has length {len(arr)}, expected {n}")
    return arr


def resolve_subscript(
    i: Any, n: int, names: np.ndarray | None
) -> np.ndarray:
    """Subscript -> int64 positions. Ints are 0-based (numpy rules apart from bounds)."""
    if isinstance(i, slice):
        return np.arange(n, dtype=np.int64)[i]
    if isinstance(i, (int, np.integer)):
        i = [i]
    arr = np.asarray(i)
    if arr.size == 0:
        return np.zeros(0, dtype=np.int64)
    if arr.dtype == bool:
        if len(arr) != n:
            raise ValueError(f"logical subscript has length {len(arr)}, expected {n}")
        return np.flatnonzero(arr).astype(np.int64)
    if arr.dtype.kind in "US" or arr.dtype == object:
        if names is None:
            raise ValueError("subscript contains invalid names")
        lookup: dict[str, int] = {}
        for pos, nm in enumerate(names.tolist()):
            if nm is not None:
                lookup.setdefault(nm, pos)
        try:
            return np.array([lookup[s] for s in arr.tolist()], dtype=np.int64)
        except KeyError:
            raise ValueError("subscript contains invalid names") from None
    if arr.dtype.kind not in "iu":
        raise TypeError(f"unsupported subscript dtype {arr.dtype}")
    arr = arr.astype(np.int64)
    if arr.min() < -n or arr.max() >= n:
        raise ValueError("subscript contains out-of-bounds indices")
    return np.where(arr < 0, arr + n, arr)


class IRanges:
    """Immutable-by-convention vector of integer ranges with names and mcols."""

    __slots__ = ("_start", "_width", "_names", "_mcols")

    def __init__(
        self,
        start: Any = None,
        end: Any = None,
        width: Any = None,
        names: Sequence[str | None] | None = None,
        mcols: Any = None,
    ) -> None:
        given = [a is not None for a in (start, end, width)]
        if sum(given) < 2:
            raise ValueError("at least two of 'start', 'end' and 'width' must be supplied")
        s, e, w = (
            as_int(a, nm) if a is not None else None
            for a, nm in ((start, "start"), (end, "end"), (width, "width"))
        )
        s, e, w = recycle([s, e, w])
        if s is not None and e is not None:
            inferred = e - s + 1
            if w is not None and (w != inferred).any():
                raise ValueError("'width' does not match width inferred from 'start' and 'end'")
            w = inferred
        elif w is not None and e is not None:
            s = e - w + 1
        assert s is not None and w is not None
        self._init(s, w, as_names(names, len(s)), _mcols.coerce(mcols, len(s)))

    def _init(self, start: np.ndarray, width: np.ndarray, names: np.ndarray | None, mcols: pl.DataFrame) -> None:
        if (width < 0).any():
            raise ValueError("each range must have a non-negative width")
        end = start + width - 1
        if start.size and (end.max() > 2**31 - 1 or start.min() < -(2**31) + 1):
            raise ValueError("range exceeds R integer bounds")
        self._start = start
        self._width = width
        self._names = names
        self._mcols = mcols

    @classmethod
    def _make(
        cls, start: np.ndarray, width: np.ndarray, names: np.ndarray | None = None, mcols: pl.DataFrame | None = None
    ) -> Self:
        obj = object.__new__(cls)
        obj._init(
            np.asarray(start, dtype=np.int64),
            np.asarray(width, dtype=np.int64),
            names,
            mcols if mcols is not None else _mcols.empty(),
        )
        return obj

    # accessors -----------------------------------------------------------
    @property
    def start(self) -> np.ndarray:
        return self._start.copy()

    @property
    def end(self) -> np.ndarray:
        return self._start + self._width - 1

    @property
    def width(self) -> np.ndarray:
        return self._width.copy()

    @property
    def names(self) -> np.ndarray | None:
        return None if self._names is None else self._names.copy()

    @property
    def mcols(self) -> pl.DataFrame:
        return self._mcols

    def with_mcols(self, mcols: Any) -> Self:
        return self._make(self._start, self._width, self._names, _mcols.coerce(mcols, len(self)))

    def with_names(self, names: Sequence[str | None] | None) -> Self:
        return self._make(self._start, self._width, as_names(names, len(self)), self._mcols)

    def __len__(self) -> int:
        return len(self._start)

    # subsetting / combining ----------------------------------------------
    def _take(self, idx: np.ndarray) -> Self:
        return self._make(
            self._start[idx],
            self._width[idx],
            None if self._names is None else self._names[idx],
            _mcols.take(self._mcols, idx),
        )

    def __getitem__(self, i: Any) -> Self:
        return self._take(resolve_subscript(i, len(self), self._names))

    @classmethod
    def concat(cls, *xs: Self) -> Self:
        if not xs:
            return cls._make(np.zeros(0, np.int64), np.zeros(0, np.int64))
        names: np.ndarray | None = None
        if any(x._names is not None for x in xs):
            # R pads unnamed pieces with "" when any piece has names
            names = np.concatenate(
                [x._names if x._names is not None else np.full(len(x), "", dtype=object) for x in xs]
            )
        return cls._make(
            np.concatenate([x._start for x in xs]),
            np.concatenate([x._width for x in xs]),
            names,
            _mcols.concat([x._mcols for x in xs], [len(x) for x in xs]),
        )

    def __repr__(self) -> str:
        rows = ", ".join(f"[{s},{e}]" for s, e in zip(self._start.tolist(), self.end.tolist(), strict=True))
        return f"IRanges({rows})"
