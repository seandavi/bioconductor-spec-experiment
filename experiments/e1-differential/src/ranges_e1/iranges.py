"""IRanges: immutable vector of 1-based, closed integer ranges."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _meta
from ._kernels import inter, intra, overlaps, setops
from .hits import Hits
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

    # -- intra-range ---------------------------------------------------------
    def _meta_kept(self, use_names: bool) -> tuple[Names | None, pl.DataFrame | None]:
        return (self._names if use_names else None), self._mcols

    def _int_arg(self, v: Any, what: str, **kw: Any) -> np.ndarray:
        """Integer-valued argument recycled to len(self); doubles truncate like as.integer."""
        a = np.asarray(v)
        if a.dtype.kind == "f" and not np.ma.isMaskedArray(v):
            if np.isnan(a).any():
                raise ValueError(f"'{what}' must be a numeric vector without NA's")
            a = np.trunc(a)
        elif a.dtype.kind not in "iufb" and not np.ma.isMaskedArray(v):
            raise ValueError(f"'{what}' must be a numeric vector")
        return _meta.recycle_arg(_meta.as_int_array(a if not np.ma.isMaskedArray(v) else v, f"'{what}'"), len(self), what, **kw)

    def shift(self, shift: Any = 0, use_names: bool = True) -> Self:
        by = self._int_arg(shift, "shift")
        s, w = intra.shift(self._start, self._width, by)
        return self._new(s, w, *self._meta_kept(use_names))

    def narrow(self, start: Any = None, end: Any = None, width: Any = None, use_names: bool = True) -> Self:
        """Narrow relative to each range; negative start/end count from the end (R)."""
        n = len(self)

        def sew(v: Any, what: str) -> np.ndarray:
            if v is None:
                return np.full(n, intra.NA, dtype=np.int64)
            if np.ma.isMaskedArray(v):
                arr = np.where(np.ma.getmaskarray(v), intra.NA, np.ma.getdata(v)).astype(np.int64)
            else:
                f = np.atleast_1d(np.asarray(v, dtype=float))
                arr = np.where(np.isnan(f), intra.NA, np.trunc(np.nan_to_num(f))).astype(np.int64)
            return _meta.recycle_arg(arr, n, what, warn=False)

        s, w = intra.narrow(self._start, self._width, sew(start, "start"), sew(end, "end"), sew(width, "width"))
        return self._new(s, w, *self._meta_kept(use_names))

    def resize(self, width: Any, fix: Any = "start", use_names: bool = True) -> Self:
        new_w = self._int_arg(width, "width", longer="truncate")
        if (new_w < 0).any():
            raise ValueError("'width' values must be non-negative")
        s, w = intra.resize(self._start, self._width, new_w, self._fix_codes(fix))
        return self._new(s, w, *self._meta_kept(use_names))

    def _fix_codes(self, fix: Any) -> np.ndarray:
        values = [fix] if isinstance(fix, str) else list(fix)
        table = {"start": intra.FIX_START, "end": intra.FIX_END, "center": intra.FIX_CENTER}
        if not values or any(v not in table for v in values):
            raise ValueError("'fix' must be a character vector with values in \"start\", \"end\", and \"center\"")
        codes = np.array([table[v] for v in values], dtype=np.int64)
        return _meta.recycle_arg(codes, len(self), "fix", longer="truncate")

    def flank(self, width: Any, start: Any = True, both: bool = False, use_names: bool = True) -> Self:
        by = self._int_arg(width, "width")
        flags = np.asarray(start)
        if flags.dtype.kind != "b":
            raise ValueError("'start' must be logical without NA's")
        at_start = _meta.recycle_arg(np.atleast_1d(flags), len(self), "start", longer="truncate")
        s, w = intra.flank(self._start, self._width, by, at_start, bool(both))
        return self._new(s, w, *self._meta_kept(use_names))

    def promoters(self, upstream: Any = 2000, downstream: Any = 200, use_names: bool = True) -> Self:
        up = self._int_arg(upstream, "upstream", longer="truncate")
        down = self._int_arg(downstream, "downstream", longer="truncate")
        s, w = intra.promoters(self._start, self._width, up, down)
        return self._new(s, w, *self._meta_kept(use_names))

    # -- inter-range ---------------------------------------------------------
    def _group0(self) -> np.ndarray:
        return np.zeros(len(self), dtype=np.int64)

    def range(self, *others: IRanges) -> Self:
        """Smallest range covering everything (zero-width ranges included)."""
        x = self.concat([self, *others]) if others else self
        red = inter.range_(x._group0(), x._start, x._width)
        return self._new(red.start, red.width)

    def reduce(
        self, drop_empty_ranges: bool = False, min_gapwidth: int = 1, with_revmap: bool = False
    ) -> Self:
        if min_gapwidth < 0:
            raise ValueError("'min.gapwidth' must be non-negative")
        red = inter.reduce(self._group0(), self._start, self._width, int(min_gapwidth), bool(drop_empty_ranges))
        mc = _meta.revmap_frame(red.members) if with_revmap else None
        return self._new(red.start, red.width, None, mc)

    def gaps(self, start: int | None = None, end: int | None = None) -> Self:
        """Uncovered stretches within [start, end]; defaults to the extent of non-empty ranges."""
        ne = self._width > 0
        have = bool(ne.any())
        lo = start if start is not None else (int(self._start[ne].min()) if have else None)
        hi = end if end is not None else (int(self.end[ne].max()) if have else None)
        if lo is None or hi is None:
            return self._new(*(np.empty(0, dtype=np.int64),) * 2)
        _, s, w = inter.gaps(
            self._group0(), self._start, self._width, np.array([lo], dtype=np.int64), np.array([hi], dtype=np.int64)
        )
        return self._new(s, w)

    def disjoin(self, with_revmap: bool = False) -> Self:
        red = inter.disjoin(self._group0(), self._start, self._width, with_revmap)
        mc = _meta.revmap_frame(red.members) if with_revmap else None
        return self._new(red.start, red.width, None, mc)

    # -- set operations ------------------------------------------------------
    def _set_op(self, other: IRanges, fn: Any) -> Self:
        _, s, w = fn(self._group0(), self._start, self._width, other._group0(), other._start, other._width)
        return self._new(s, w)

    def union(self, other: IRanges) -> Self:
        return self._set_op(other, setops.union)

    def intersect(self, other: IRanges) -> Self:
        return self._set_op(other, setops.intersect)

    def setdiff(self, other: IRanges) -> Self:
        return self._set_op(other, setops.setdiff)

    def _check_parallel(self, other: IRanges) -> None:
        if len(self) != len(other):
            raise ValueError("'x' and 'y' must have the same length")

    def _parallel_result(self, s: np.ndarray, w: np.ndarray) -> Self:
        """p* results keep x's names; mcols are dropped."""
        return self._new(s, w, self._names, None)

    def punion(self, other: IRanges, fill_gap: bool = False) -> Self:
        self._check_parallel(other)
        if not isinstance(fill_gap, (bool, np.bool_)):
            raise ValueError("'fill.gap' must be TRUE or FALSE")
        return self._parallel_result(*setops.punion(self._start, self._width, other._start, other._width, bool(fill_gap)))

    def pintersect(self, other: IRanges, resolve_empty: str = "none") -> Self:
        self._check_parallel(other)
        return self._parallel_result(
            *setops.pintersect(self._start, self._width, other._start, other._width, resolve_empty)
        )

    def psetdiff(self, other: IRanges) -> Self:
        self._check_parallel(other)
        return self._parallel_result(*setops.psetdiff(self._start, self._width, other._start, other._width))

    # -- overlaps ------------------------------------------------------------
    def _pairs(self, subject: IRanges, maxgap: int, minoverlap: int, type_: str) -> tuple[np.ndarray, np.ndarray]:
        """Hit pairs for (self as query, subject); subclasses add compatibility filters."""
        return overlaps.find(self._start, self._width, subject._start, subject._width, maxgap, minoverlap, type_)

    def find_overlaps(
        self,
        subject: IRanges | None = None,
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",  # noqa: A002 - R's argument name
        select: str = "all",
        drop_self: bool = False,
        drop_redundant: bool = False,
    ) -> Hits | np.ma.MaskedArray:
        """Overlaps of self (query) with `subject`; `subject=None` means self-overlaps.

        select="all" -> Hits; "first" / "last" / "arbitrary" -> masked array of 0-based
        subject indices (masked = no hit).
        """
        overlaps.check_args(type, int(maxgap), int(minoverlap))
        if select not in ("all", "first", "last", "arbitrary"):
            raise ValueError("'arg' should be one of \"all\", \"first\", \"last\", \"arbitrary\"")
        self_hits = subject is None
        subj = self if subject is None else subject
        q, s = self._pairs(subj, int(maxgap), int(minoverlap), type)
        if self_hits:
            keep = np.ones(len(q), dtype=bool)
            if drop_self:
                keep &= q != s
            if drop_redundant:  # drop (i, j), i > j, when its mirror (j, i) is also a hit
                n = len(self)
                key = set((q * n + s).tolist())
                mirror = np.array([(b * n + a) in key for a, b in zip(q.tolist(), s.tolist())], dtype=bool)
                keep &= ~((q > s) & mirror)
            q, s = q[keep], s[keep]
        hits = Hits(q, s, len(self), len(subj))
        if select == "arbitrary" and overlaps.uses_query_nclist(len(self), len(subj)):
            select = "last"  # R's engine reports the last-seen subject in this mode
        return hits if select == "all" else hits.select(select)

    def count_overlaps(
        self, subject: IRanges | None = None, maxgap: int = -1, minoverlap: int = 0, type: str = "any"  # noqa: A002
    ) -> np.ndarray:
        hits = self.find_overlaps(subject, maxgap, minoverlap, type)
        assert isinstance(hits, Hits)
        return hits.count_by_query()

    def overlaps_any(
        self, subject: IRanges | None = None, maxgap: int = -1, minoverlap: int = 0, type: str = "any"  # noqa: A002
    ) -> np.ndarray:
        return self.count_overlaps(subject, maxgap, minoverlap, type) > 0

    def subset_by_overlaps(
        self,
        subject: IRanges,
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",  # noqa: A002
        invert: bool = False,
    ) -> Self:
        hit = self.overlaps_any(subject, maxgap, minoverlap, type)
        return self[~hit if invert else hit]
