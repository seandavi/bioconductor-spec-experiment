"""IRanges: immutable vector of 1-based, closed integer ranges."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _meta
from ._meta import INT_MAX, Names


def _solve_sew(
    start: Any, end: Any, width: Any
) -> tuple[np.ndarray, np.ndarray]:
    """R's `solveUserSEW`: any two of start/end/width -> (start, width)."""
    given = [a is not None for a in (start, end, width)]
    if sum(given) < 2:
        raise ValueError("at least two of the 'start', 'end', and 'width' arguments must be supplied")
    s = _meta.as_int_array(start, "'start'") if start is not None else None
    e = _meta.as_int_array(end, "'end'") if end is not None else None
    w = _meta.as_int_array(width, "'width'") if width is not None else None
    s, e, w = _meta.recycle([s, e, w])
    if s is not None and e is not None:
        inferred = e - s + 1
        if w is not None and (w != inferred).any():
            k = int(np.flatnonzero(w != inferred)[0])
            raise ValueError(
                f"In range {k + 1}: the supplied 'width' ({w[k]}) doesn't match the width "
                f"inferred from the supplied 'start' and 'end' ({inferred[k]})."
            )
        w = inferred
    elif s is not None and w is not None:
        pass
    else:
        assert e is not None and w is not None
        s = e - w + 1
    assert s is not None and w is not None
    return s, w


def _runs(x: Any) -> tuple[np.ndarray, np.ndarray]:
    """`IRanges(start=<vector>)`: logical -> runs of TRUE; integer -> runs of +1 steps."""
    a = np.asarray(x)
    if a.dtype.kind == "b":
        idx = np.flatnonzero(a) + 1
    else:
        idx = _meta.as_int_array(x, "'start'")
    if len(idx) == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    breaks = np.flatnonzero(np.diff(idx) != 1) + 1
    first = np.concatenate([[0], breaks])
    last = np.concatenate([breaks, [len(idx)]]) - 1
    return idx[first], idx[last] - idx[first] + 1


def _validate(start: np.ndarray, width: np.ndarray) -> None:
    if (width < 0).any():
        raise ValueError("each range must have a non-negative width")
    end = start + width - 1
    if len(start) and (
        (start > INT_MAX).any() or (start < -INT_MAX).any() or (end > INT_MAX).any() or (end < -INT_MAX).any()
    ):
        raise ValueError("each range must have an end that is < 2^31 and > - 2^31")


class IRanges:
    """Integer ranges; coordinates 1-based closed, `end = start + width - 1`."""

    __slots__ = ("_start", "_width", "_names", "_mcols")

    def __init__(
        self,
        start: Any = None,
        end: Any = None,
        width: Any = None,
        names: Iterable[str | None] | None = None,
        mcols: Any = None,
    ) -> None:
        if start is not None and end is None and width is None:
            s, w = _runs(start)
        else:
            s, w = _solve_sew(start, end, width)
        _validate(s, w)
        self._init(s, w, _meta.check_names(names, len(s)), _meta.check_mcols(mcols, len(s)))

    def _init(
        self, start: np.ndarray, width: np.ndarray, names: Names | None, mcols: pl.DataFrame | None
    ) -> None:
        start = np.ascontiguousarray(start, dtype=np.int64)
        width = np.ascontiguousarray(width, dtype=np.int64)
        start.flags.writeable = False
        width.flags.writeable = False
        self._start = start
        self._width = width
        self._names = names
        self._mcols = mcols

    @classmethod
    def _new(
        cls,
        start: np.ndarray,
        width: np.ndarray,
        names: Names | None = None,
        mcols: pl.DataFrame | None = None,
        *,
        validate: bool = True,
    ) -> Self:
        """Internal constructor from (start, width) arrays."""
        start = np.asarray(start, dtype=np.int64)
        width = np.asarray(width, dtype=np.int64)
        if validate:
            _validate(start, width)
        obj = cls.__new__(cls)
        obj._init(start, width, names, mcols)
        return obj

    # -- accessors ---------------------------------------------------------
    @property
    def start(self) -> np.ndarray:
        return self._start

    @property
    def width(self) -> np.ndarray:
        return self._width

    @property
    def end(self) -> np.ndarray:
        e = self._start + self._width - 1
        e.flags.writeable = False
        return e

    @property
    def names(self) -> Names | None:
        return self._names

    @property
    def mcols(self) -> pl.DataFrame | None:
        return self._mcols

    def __len__(self) -> int:
        return len(self._start)

    def __repr__(self) -> str:
        n = len(self)
        rows = ", ".join(f"{s}-{s + w - 1}" for s, w in zip(self._start[:6], self._width[:6]))
        more = ", ..." if n > 6 else ""
        return f"IRanges(n={n}: {rows}{more})"

    def __iter__(self) -> Iterable[tuple[int, int]]:
        return iter(zip(self._start.tolist(), self.end.tolist()))

    # -- setters (copies) --------------------------------------------------
    def set_names(self, names: Iterable[str | None] | None) -> Self:
        return self._new(self._start, self._width, _meta.check_names(names, len(self)), self._mcols, validate=False)

    def set_mcols(self, mcols: Any) -> Self:
        return self._new(self._start, self._width, self._names, _meta.check_mcols(mcols, len(self)), validate=False)

    def _with(self, start: np.ndarray, width: np.ndarray, *, keep_meta: bool = True) -> Self:
        """New object with new coordinates; names/mcols kept (same length) or dropped."""
        return self._new(
            start, width, self._names if keep_meta else None, self._mcols if keep_meta else None
        )

    # -- subsetting / combining --------------------------------------------
    def _take(self, pos: np.ndarray) -> Self:
        return self._new(
            self._start[pos],
            self._width[pos],
            _meta.take_names(self._names, pos),
            _meta.take_mcols(self._mcols, pos),
            validate=False,
        )

    def __getitem__(self, i: Any) -> Self:
        return self._take(_meta.resolve_index(i, len(self), self._names))

    @classmethod
    def concat(cls, parts: Sequence[IRanges]) -> Self:
        parts = list(parts)
        if not parts or not all(isinstance(p, IRanges) for p in parts):
            raise TypeError("concat needs a non-empty sequence of IRanges")
        return cls._new(
            np.concatenate([p._start for p in parts]),
            np.concatenate([p._width for p in parts]),
            _meta.concat_names([(p._names, len(p)) for p in parts]),
            _meta.concat_mcols([(p._mcols, len(p)) for p in parts]),
        )
