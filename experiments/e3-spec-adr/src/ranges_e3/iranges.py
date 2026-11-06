"""``IRanges``: integer ranges, 1-based and closed (spec 10-iranges-core.md)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from ranges_e3 import _kernels, _mcols, _subscript
from ranges_e3._types import BoolArray, IntArray
from ranges_e3.hits import Hits
from ranges_e3.rle import Rle


def as_int_array(values: Any, what: str) -> IntArray:
    """1-D int64 array from ints/sequences; NA and non-integers refused (CORE.4)."""
    if values is None:
        raise ValueError(f"{what} is missing")
    if isinstance(values, np.ma.MaskedArray) and np.ma.is_masked(values):
        raise ValueError(f"{what} contains NA")
    arr = np.atleast_1d(np.asarray(values))
    if arr.ndim != 1:
        raise ValueError(f"{what} must be 1-dimensional")
    if arr.dtype == object:
        if any(v is None for v in arr):
            raise ValueError(f"{what} contains NA")
        arr = arr.astype(np.float64)
    if arr.dtype.kind == "f":
        if np.isnan(arr).any():
            raise ValueError(f"{what} contains NA")
        if (arr != np.floor(arr)).any():
            raise ValueError(f"{what} must be integers")
    elif arr.dtype.kind not in "iub":
        raise ValueError(f"{what} must be integers")
    return arr.astype(np.int64)


def per_range(values: Any, n: int, what: str) -> IntArray:
    """Broadcast a scalar / length-1 / length-``n`` argument to length ``n`` (CONV.6)."""
    arr = as_int_array(values, what)
    if arr.size == 1:
        return np.full(n, arr[0], dtype=np.int64)
    if arr.size != n:
        raise ValueError(f"{what} has length {arr.size}, expected 1 or {n}")
    return arr


class IRanges:
    """Immutable vector of closed integer ranges with optional names and mcols."""

    __slots__ = ("_start", "_width", "_names", "_mcols")

    def __init__(
        self,
        start: Any = None,
        end: Any = None,
        width: Any = None,
        names: Sequence[str | None] | None = None,
        mcols: pl.DataFrame | None = None,
    ) -> None:
        given = [v is not None for v in (start, end, width)]
        if sum(given) != 2:
            raise ValueError("exactly two of start, end, width must be given")
        parts = {
            k: as_int_array(v, k)
            for k, v in (("start", start), ("end", end), ("width", width))
            if v is not None
        }
        n = max(len(a) for a in parts.values())
        if any(len(a) not in (1, n) for a in parts.values()):
            raise ValueError("start/end/width lengths do not match")
        parts = {k: np.broadcast_to(a, n).copy() for k, a in parts.items()}
        if start is None:
            s = parts["end"] - parts["width"] + 1
            w = parts["width"]
        elif end is None:
            s, w = parts["start"], parts["width"]
        else:
            s = parts["start"]
            w = parts["end"] - s + 1
        if (w < 0).any():
            raise ValueError("negative width: each range needs end >= start - 1")
        self._start: IntArray = s
        self._width: IntArray = w
        self._names = self._check_names(names, n)
        self._mcols = _mcols.validate(mcols, n)

    @staticmethod
    def _check_names(names: Sequence[str | None] | None, n: int) -> list[str | None] | None:
        if names is None:
            return None
        out = list(names)
        if len(out) != n:
            raise ValueError(f"names has length {len(out)}, expected {n}")
        return out

    @classmethod
    def _new(
        cls,
        start: IntArray,
        width: IntArray,
        names: list[str | None] | None = None,
        mcols: pl.DataFrame | None = None,
    ) -> Self:
        """Trusted constructor: arrays are already validated int64, width >= 0."""
        obj = cls.__new__(cls)
        obj._start = start
        obj._width = width
        obj._names = names
        obj._mcols = mcols
        return obj

    # --- accessors (CORE.6-8) ---

    def __len__(self) -> int:
        return len(self._start)

    @property
    def start(self) -> IntArray:
        return self._start.copy()

    @property
    def width(self) -> IntArray:
        return self._width.copy()

    @property
    def end(self) -> IntArray:
        return self._start + self._width - 1

    @property
    def names(self) -> list[str | None] | None:
        return None if self._names is None else list(self._names)

    @property
    def mcols(self) -> pl.DataFrame:
        return _mcols.public(self._mcols)

    def with_names(self, names: Sequence[str | None] | None) -> Self:
        return self._new(self._start, self._width, self._check_names(names, len(self)), self._mcols)

    def with_mcols(self, mcols: pl.DataFrame | None) -> Self:
        return self._new(self._start, self._width, self._names, _mcols.validate(mcols, len(self)))

    def bare(self) -> Self:
        """Same ranges without names and mcols."""
        return self._new(self._start, self._width)

    # --- subsetting (CORE.9-10) ---

    def __getitem__(self, index: Any) -> Self:
        pos = _subscript.resolve(index, len(self), self._names)
        return self._take(pos)

    def _take(self, pos: IntArray) -> Self:
        names = None if self._names is None else [self._names[i] for i in pos]
        return self._new(self._start[pos], self._width[pos], names, _mcols.take(self._mcols, pos))

    # --- intra-range (spec 20-intra.md) ---

    def _like(self, start: IntArray, width: IntArray) -> Self:
        """Same names/mcols, new coordinates."""
        return self._new(start, width, self._names, self._mcols)

    def shift(self, shift: Any) -> Self:
        """INTRA.1."""
        return self._like(*_kernels.shift(self._start, self._width, per_range(shift, len(self), "shift")))

    def narrow(self, start: Any = None, end: Any = None, width: Any = None) -> Self:
        """INTRA.2 (positions relative to each range)."""
        n = len(self)
        args = [None if v is None else per_range(v, n, k) for k, v in (("start", start), ("end", end), ("width", width))]
        return self._like(*_kernels.narrow(self._start, self._width, *args))

    def resize(self, width: Any, fix: str = "start") -> Self:
        """INTRA.3."""
        code = {"start": _kernels.FIX_START, "end": _kernels.FIX_END, "center": _kernels.FIX_CENTER}
        if fix not in code:
            raise ValueError("fix must be 'start', 'end' or 'center'")
        n = len(self)
        fixes = np.full(n, code[fix], dtype=np.int64)
        return self._like(*_kernels.resize(self._start, self._width, per_range(width, n, "width"), fixes))

    def flank(self, width: Any, start: Any = True, both: bool = False) -> Self:
        """INTRA.4."""
        n = len(self)
        at_start = np.broadcast_to(np.asarray(start, dtype=np.bool_), n)
        return self._like(*_kernels.flank(self._start, self._width, per_range(width, n, "width"), at_start, both))

    def promoters(self, upstream: Any = 2000, downstream: Any = 200) -> Self:
        """INTRA.5."""
        n = len(self)
        mirror = np.zeros(n, dtype=np.bool_)
        return self._like(
            *_kernels.promoters(
                self._start, self._width, per_range(upstream, n, "upstream"), per_range(downstream, n, "downstream"), mirror
            )
        )

    # --- inter-range (spec 30-inter.md) ---

    @property
    def _end(self) -> IntArray:
        return self._start + self._width - 1

    def range(self) -> Self:
        """INTER.1."""
        bounds = _kernels.span(self._start, self._end)
        if bounds is None:
            return self._new(_mcols.empty_int(), _mcols.empty_int())
        return self._new(np.array([bounds[0]], dtype=np.int64), np.array([bounds[1] - bounds[0] + 1], dtype=np.int64))

    def reduce(self, min_gapwidth: int = 1, drop_empty_ranges: bool = False, with_revmap: bool = False) -> Self:
        """INTER.2-5."""
        gs, ge, order, bounds = _kernels.reduce_scan(self._start, self._end, int(min_gapwidth), drop_empty_ranges)
        out = self._new(gs, ge - gs + 1)
        if with_revmap:
            out = out.with_mcols(pl.DataFrame([_mcols.revmap_series(_kernels.split_groups(order, bounds))]))
        return out

    def gaps(self, start: int | None = None, end: int | None = None) -> Self:
        """INTER.6."""
        bounds = _kernels.span(self._start, self._end)
        lo = start if start is not None else (None if bounds is None else bounds[0])
        hi = end if end is not None else (None if bounds is None else bounds[1])
        if lo is None or hi is None:
            return self._new(_mcols.empty_int(), _mcols.empty_int())
        gs, ge = _kernels.gaps(self._start, self._end, int(lo), int(hi))
        return self._new(gs, ge - gs + 1)

    def disjoin(self, with_revmap: bool = False) -> Self:
        """INTER.7-8."""
        ps, pe, revmap = _kernels.disjoin(self._start, self._end, with_revmap)
        out = self._new(ps, pe - ps + 1)
        if with_revmap:
            out = out.with_mcols(pl.DataFrame([_mcols.revmap_series(revmap)]))
        return out

    # --- set operations (spec 40-setops.md) ---

    def _set_op(self, y: IRanges, op: str) -> Self:
        s, e = _kernels.set_operation(self._start, self._end, y._start, y._end, op)
        return self._new(s, e - s + 1)

    def union(self, y: IRanges) -> Self:
        """SET.1."""
        return self._set_op(y, "union")

    def intersect(self, y: IRanges) -> Self:
        """SET.1."""
        return self._set_op(y, "intersect")

    def setdiff(self, y: IRanges) -> Self:
        """SET.1."""
        return self._set_op(y, "setdiff")

    def _parallel(self, y: IRanges) -> tuple[IntArray, IntArray, IntArray, IntArray]:
        if len(self) != len(y):
            raise ValueError("parallel operations need equal lengths (SET.2)")
        return self._start, self._end, y._start, y._end

    def _parallel_result(self, s: IntArray, e: IntArray) -> Self:
        return self._new(s, e - s + 1, self._names)

    def punion(self, y: IRanges, fill_gap: bool = False) -> Self:
        """SET.3."""
        return self._parallel_result(*_kernels.punion(*self._parallel(y), fill_gap))

    def pintersect(self, y: IRanges, resolve_empty: str = "none") -> Self:
        """SET.4."""
        return self._parallel_result(*_kernels.pintersect(*self._parallel(y), resolve_empty))

    def psetdiff(self, y: IRanges) -> Self:
        """SET.5."""
        return self._parallel_result(*_kernels.psetdiff(*self._parallel(y)))

    # --- coverage (spec 50-coverage.md) ---

    def coverage(self, shift: Any = 0, width: int | None = None) -> Rle:
        """COV.1-3."""
        if width is not None and width < 0:
            raise ValueError("width must be non-negative")
        lengths, values = _kernels.coverage_runs(
            self._start, self._end, per_range(shift, len(self), "shift"), width
        )
        return Rle(lengths, values)

    # --- overlaps (spec 60-overlaps.md) ---

    def _overlap_pairs(
        self, subject: IRanges, max_gap: int, min_overlap: int, type: str
    ) -> tuple[IntArray, IntArray]:
        return _kernels.find_overlaps(
            self._start, self._end, subject._start, subject._end, int(max_gap), int(min_overlap), type
        )

    def find_overlaps(
        self, subject: IRanges, max_gap: int = -1, min_overlap: int = 0, type: str = "any", select: str = "all"
    ) -> Hits | np.ma.MaskedArray:  # type: ignore[type-arg]
        """OVL.2-9."""
        qi, si = self._overlap_pairs(subject, max_gap, min_overlap, type)
        if select == "all":
            return Hits(qi, si, len(self), len(subject))
        return _kernels.select_pairs(qi, si, len(self), select)

    def count_overlaps(self, subject: IRanges, max_gap: int = -1, min_overlap: int = 0, type: str = "any") -> IntArray:
        """OVL.10."""
        qi, _ = self._overlap_pairs(subject, max_gap, min_overlap, type)
        return np.bincount(qi, minlength=len(self)).astype(np.int64)

    def overlaps_any(self, subject: IRanges, max_gap: int = -1, min_overlap: int = 0, type: str = "any") -> BoolArray:
        """OVL.10."""
        return self.count_overlaps(subject, max_gap, min_overlap, type) > 0

    def subset_by_overlaps(
        self, subject: IRanges, max_gap: int = -1, min_overlap: int = 0, type: str = "any", invert: bool = False
    ) -> Self:
        """OVL.10."""
        hit = self.overlaps_any(subject, max_gap, min_overlap, type)
        return self._take(np.flatnonzero(~hit if invert else hit).astype(np.int64))

    def __repr__(self) -> str:
        shown = ", ".join(f"[{s},{e}]" for s, e in zip(self._start[:6], self.end[:6]))
        more = ", ..." if len(self) > 6 else ""
        return f"IRanges(n={len(self)}: {shown}{more})"


def concat_iranges(objs: Sequence[IRanges]) -> IRanges:
    """Append ranges in order (CORE.11-12)."""
    if not objs:
        raise ValueError("concat needs at least one object")
    start = np.concatenate([o._start for o in objs])
    width = np.concatenate([o._width for o in objs])
    names: list[str | None] | None = None
    if any(o._names is not None for o in objs):
        names = []
        for o in objs:
            names.extend(o._names if o._names is not None else [""] * len(o))
    mcols = _mcols.concat([(o._mcols, len(o)) for o in objs])
    return IRanges._new(start, width, names, mcols)
