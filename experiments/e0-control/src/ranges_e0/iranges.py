"""IRanges: integer ranges, 1-based and closed (``end = start + width - 1``)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import polars as pl

from . import _intra, _order, _validate as V

NO_MATCH: int = -1
"""Marker for "no match" in index-valued results (0-based indices otherwise).

Beware: ``-1`` is also a valid Python index (the last element)."""


class IRanges:
    """A vector of integer ranges with optional names and metadata columns.

    Immutable by convention: every operation returns a new object. Index-valued
    results are 0-based (see ``NO_MATCH``). Mirrors R's ``IRanges``.
    """

    __slots__ = ("_start", "_width", "_names", "_mcols")

    def __init__(
        self,
        start: Any = None,
        end: Any = None,
        width: Any = None,
        names: Sequence[str] | None = None,
        mcols: pl.DataFrame | dict[str, Sequence[Any]] | None = None,
    ) -> None:
        s = None if start is None else V.as_int_array(start, "start")
        e = None if end is None else V.as_int_array(end, "end")
        w = None if width is None else V.as_int_array(width, "width")
        given = [a for a in (s, e, w) if a is not None]
        if s is not None and e is None and w is None:
            w = np.ones(1, dtype=np.int64)
            given = [s, w]
        if len(given) < 2:
            raise ValueError("at least two of 'start', 'end', and 'width' must be supplied")
        n = max(len(a) for a in given)
        if any(len(a) == 0 for a in given):
            n = 0
        s = None if s is None else V.recycle(s, n, "start")
        e = None if e is None else V.recycle(e, n, "end")
        w = None if w is None else V.recycle(w, n, "width")
        if s is not None and e is not None:
            inferred = e - s + 1
            if w is not None and np.any(w != inferred):
                raise ValueError(
                    "the supplied 'width' doesn't match the width inferred from 'start' and 'end'"
                )
            w = inferred
        elif s is not None:
            assert w is not None
        else:
            assert e is not None and w is not None
            s = e - w + 1
        assert s is not None and w is not None
        if np.any(w < 0):
            raise ValueError("each range must have an end that is greater or equal to its start minus one")
        self._start = s
        self._width = w
        self._names = V.as_names(names, n)
        self._mcols = V.as_mcols(mcols, n)

    @classmethod
    def _make(
        cls,
        start: np.ndarray,
        width: np.ndarray,
        names: np.ndarray | None = None,
        mcols: pl.DataFrame | None = None,
    ) -> IRanges:
        """Internal constructor from validated int64 arrays (no checks except width)."""
        if len(width) and width.min() < 0:
            raise ValueError("negative width")
        obj = object.__new__(cls)
        obj._start = start
        obj._width = width
        obj._names = names
        obj._mcols = pl.DataFrame() if mcols is None else mcols
        return obj

    @classmethod
    def empty(cls) -> IRanges:
        return cls._make(np.array([], dtype=np.int64), np.array([], dtype=np.int64))

    # ----- accessors -------------------------------------------------------
    @property
    def start(self) -> np.ndarray:
        return _ro(self._start)

    @property
    def width(self) -> np.ndarray:
        return _ro(self._width)

    @property
    def end(self) -> np.ndarray:
        return self._start + self._width - 1

    @property
    def names(self) -> list[str] | None:
        return None if self._names is None else list(self._names)

    @property
    def mcols(self) -> pl.DataFrame:
        return self._mcols

    def __len__(self) -> int:
        return len(self._start)

    def with_names(self, names: Sequence[str] | None) -> IRanges:
        """R's ``names<-``; ``None`` removes names."""
        return IRanges._make(self._start, self._width, V.as_names(names, len(self)), self._mcols)

    def with_mcols(self, mcols: pl.DataFrame | dict[str, Sequence[Any]] | None) -> IRanges:
        """R's ``mcols<-``; ``None`` removes all metadata columns."""
        return IRanges._make(self._start, self._width, self._names, V.as_mcols(mcols, len(self)))

    def _like(self, start: np.ndarray, width: np.ndarray, keep_names: bool = True) -> IRanges:
        """New IRanges, same length, keeping names (optional) and mcols."""
        return IRanges._make(start, width, self._names if keep_names else None, self._mcols)

    # ----- intra-range ------------------------------------------------------
    def shift(self, shift: Any, use_names: bool = True) -> IRanges:
        """R ``shift``: move every range by ``shift`` (scalar or per-range)."""
        s, w = _intra.shift_k(self._start, self._width, shift)
        return self._like(s, w, use_names)

    def narrow(
        self, start: Any = None, end: Any = None, width: Any = None, use_names: bool = True
    ) -> IRanges:
        """R ``narrow``: ``start``/``end`` are positions relative to each range
        (1 = first base; negative counts from the end). Supply ``start``
        and/or ``end``, or exactly one of them with ``width``."""
        s, w = _intra.narrow_k(self._start, self._width, start, end, width)
        return self._like(s, w, use_names)

    def resize(self, width: Any, fix: str | Sequence[str] = "start", use_names: bool = True) -> IRanges:
        """R ``resize``: set widths, anchoring at ``fix`` = start / end / center."""
        s, w = _intra.resize_k(self._start, self._width, width, _intra.as_fix(fix, len(self)))
        return self._like(s, w, use_names)

    def flank(self, width: Any, start: bool | Sequence[bool] = True, both: bool = False, use_names: bool = True) -> IRanges:
        """R ``flank``: ranges adjacent to the start (or end) of each range.
        Negative ``width`` flanks inside the range."""
        at = np.asarray(start, dtype=bool)
        at = np.repeat(at, len(self)) if at.ndim == 0 or len(at) == 1 else at
        s, w = _intra.flank_k(self._start, self._width, width, at, both)
        return self._like(s, w, use_names)

    def promoters(self, upstream: Any = 2000, downstream: Any = 200, use_names: bool = True) -> IRanges:
        """R ``promoters``: ``[start - upstream, start + downstream - 1]``."""
        s, w = _intra.promoters_k(self._start, self._width, upstream, downstream)
        return self._like(s, w, use_names)

    # ----- ordering ---------------------------------------------------------
    def _order_keys(self) -> list[np.ndarray]:
        return [self._start, self._width]

    def order(self, decreasing: bool = False) -> np.ndarray:
        """R ``order``: 0-based permutation sorting by start, then width."""
        return _order.order_k(self._order_keys(), decreasing)

    def sort(self, decreasing: bool = False) -> IRanges:
        return self._take(self.order(decreasing))

    def rank(self, ties_method: str = "first") -> np.ndarray:
        """R ``rank``: 1-based ranks. Deviation: R requires an explicit
        ``ties.method`` for IRanges; here the default is ``"first"``."""
        return _order.rank_k(self._order_keys(), ties_method)

    def duplicated(self, from_last: bool = False) -> np.ndarray:
        """R ``duplicated``: equality on (start, width); names/mcols ignored."""
        return _order.duplicated_k(self._order_keys(), from_last)

    def unique(self, from_last: bool = False) -> IRanges:
        return self._take(np.flatnonzero(~self.duplicated(from_last)))

    # ----- subsetting / combining ----------------------------------------
    def __getitem__(self, index: Any) -> IRanges:
        idx = V.resolve_index(index, len(self), self._names)
        return self._take(idx)

    def _take(self, idx: np.ndarray) -> IRanges:
        return IRanges._make(
            self._start[idx],
            self._width[idx],
            None if self._names is None else self._names[idx],
            V.take_rows(self._mcols, idx),
        )

    @classmethod
    def concat(cls, *objs: IRanges) -> IRanges:
        """R's ``c()``: names are padded with ``""``; mcols aligned by column name."""
        return concat(*objs)

    # ----- comparison / display ------------------------------------------
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, IRanges):
            return NotImplemented
        if not (np.array_equal(self._start, other._start) and np.array_equal(self._width, other._width)):
            return False
        if (self._names is None) != (other._names is None):
            return False
        if self._names is not None and other._names is not None and list(self._names) != list(other._names):
            return False
        return self._mcols.equals(other._mcols) if (self._mcols.width or other._mcols.width) else True

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        n = len(self)
        head = f"IRanges(n={n}"
        if self._mcols.width:
            head += f", mcols={self._mcols.columns}"
        rows = []
        for i in range(min(n, 10)):
            nm = f"{self._names[i]}: " if self._names is not None else ""
            rows.append(f"{nm}{self._start[i]}-{self._start[i] + self._width[i] - 1} (w={self._width[i]})")
        more = ", ..." if n > 10 else ""
        return head + (f", [{'; '.join(rows)}{more}])" if rows else ")")

    def to_tuples(self) -> list[tuple[int, int]]:
        """``[(start, end), ...]`` as plain ints (handy for tests)."""
        return list(zip(self._start.tolist(), self.end.tolist()))


def _ro(a: np.ndarray) -> np.ndarray:
    v = a.view()
    v.flags.writeable = False
    return v


def concat(*objs: Any) -> Any:
    """R's ``c()`` for IRanges or GRanges objects (dispatches on the first)."""
    if not objs:
        raise ValueError("nothing to concatenate")
    first = objs[0]
    if not isinstance(first, IRanges):
        return first.concat(*objs)
    if not all(isinstance(o, IRanges) for o in objs):
        raise TypeError("can only concatenate IRanges with IRanges")
    start = np.concatenate([o._start for o in objs])
    width = np.concatenate([o._width for o in objs])
    names = None
    if any(o._names is not None for o in objs):
        names = np.concatenate(
            [o._names if o._names is not None else np.full(len(o), "", dtype=object) for o in objs]
        )
    mcols = V.concat_mcols([o._mcols for o in objs], [len(o) for o in objs])
    return IRanges._make(start, width, names, mcols)
