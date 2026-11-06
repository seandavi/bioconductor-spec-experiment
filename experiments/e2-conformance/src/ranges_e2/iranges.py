"""IRanges: 1-based closed integer ranges with optional names and metadata columns."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import polars as pl

from . import _coverage, _inter, _intra, _mcols, _nearest, _ordering, _overlaps
from ._util import IntArray, as_bool_array, as_int64, check_int32, positions_from_subscript, rep_len
from .hits import Hits
from .rle import Rle

Names = list[str | None]


def _int_arg(v: Any, what: str) -> IntArray | None:
    return None if v is None else as_int64(v, what)


def _solve_start_width(
    start: IntArray | None, end: IntArray | None, width: IntArray | None
) -> tuple[IntArray, IntArray]:
    given = [a for a in (start, end, width) if a is not None]
    if not given:
        z = np.zeros(0, dtype=np.int64)
        return z, z
    if len(given) == 1:
        raise ValueError("at least two of start, end and width must be supplied")
    n = max(len(a) for a in given)
    if min(len(a) for a in given) == 0 and n > 0:
        raise ValueError("cannot recycle a zero-length argument")
    start, end, width = (None if a is None else rep_len(a, n) for a in (start, end, width))
    if start is not None and width is not None:
        if end is not None and np.any(end - start + 1 != width):
            raise ValueError("the supplied width does not match the width inferred from start and end")
    elif start is not None and end is not None:
        width = end - start + 1
    else:
        assert end is not None and width is not None
        start = end - width + 1
    assert start is not None and width is not None
    if np.any(width < 0):
        raise ValueError("each range must have a non-negative width (end >= start - 1)")
    return start, width


class IRanges:
    """Immutable vector of integer ranges. Operations return new objects.

    Coordinates are 1-based and closed: ``end = start + width - 1``; zero-width ranges have
    ``end = start - 1``. Indices taken and returned by methods are 0-based.
    """

    __slots__ = ("_start", "_width", "_names", "_mcols")
    _start: npt.NDArray[np.int32]
    _width: npt.NDArray[np.int32]
    _names: Names | None
    _mcols: pl.DataFrame | None

    def __init__(
        self,
        start: Any = None,
        end: Any = None,
        width: Any = None,
        names: Sequence[str | None] | None = None,
        mcols: pl.DataFrame | None = None,
    ) -> None:
        s, w = _solve_start_width(_int_arg(start, "start"), _int_arg(end, "end"), _int_arg(width, "width"))
        object.__setattr__(self, "_start", check_int32(s, "start").astype(np.int32))
        object.__setattr__(self, "_width", check_int32(w, "width").astype(np.int32))
        object.__setattr__(self, "_names", None)
        object.__setattr__(self, "_mcols", None)
        if names is not None:
            object.__setattr__(self, "_names", self._pad_names(names, len(s)))
        if mcols is not None:
            object.__setattr__(self, "_mcols", _mcols.normalise(mcols, len(s)))

    def __setattr__(self, key: str, value: Any) -> None:
        raise AttributeError("IRanges is immutable")

    @staticmethod
    def _pad_names(names: Sequence[str | None], n: int) -> Names:
        lst = list(names)
        if len(lst) > n:
            raise ValueError(f"attempt to set too many names ({len(lst)}) on IRanges object of length {n}")
        return lst + [None] * (n - len(lst))

    @classmethod
    def _make(
        cls, start: IntArray, width: IntArray, names: Names | None = None, mcols: pl.DataFrame | None = None
    ) -> "IRanges":
        """Trusted constructor from validated int64 arrays."""
        obj = object.__new__(cls)
        object.__setattr__(obj, "_start", check_int32(np.asarray(start, dtype=np.int64), "start").astype(np.int32))
        object.__setattr__(obj, "_width", check_int32(np.asarray(width, dtype=np.int64), "width").astype(np.int32))
        object.__setattr__(obj, "_names", names)
        object.__setattr__(obj, "_mcols", mcols)
        return obj

    # ---- accessors -------------------------------------------------------

    def __len__(self) -> int:
        return len(self._start)

    @property
    def start(self) -> npt.NDArray[np.int32]:
        return self._start

    @property
    def width(self) -> npt.NDArray[np.int32]:
        return self._width

    @property
    def end(self) -> npt.NDArray[np.int32]:
        return (self._start.astype(np.int64) + self._width - 1).astype(np.int32)

    @property
    def names(self) -> Names | None:
        return None if self._names is None else list(self._names)

    @property
    def mcols(self) -> pl.DataFrame:
        """Metadata columns; an object without columns returns an empty frame (height 0)."""
        return _mcols.empty_frame() if self._mcols is None else self._mcols

    @property
    def has_mcols(self) -> bool:
        return self._mcols is not None

    def with_names(self, names: Sequence[str | None] | None) -> "IRanges":
        padded = None if names is None else self._pad_names(names, len(self))
        return IRanges._make(self._start.astype(np.int64), self._width.astype(np.int64), padded, self._mcols)

    def with_mcols(self, mcols: pl.DataFrame | None) -> "IRanges":
        return IRanges._make(
            self._start.astype(np.int64),
            self._width.astype(np.int64),
            self._names,
            _mcols.normalise(mcols, len(self)),
        )

    def _s(self) -> IntArray:
        return self._start.astype(np.int64)

    def _w(self) -> IntArray:
        return self._width.astype(np.int64)

    def _e(self) -> IntArray:
        return self._s() + self._w() - 1

    # ---- subsetting / combining -----------------------------------------

    def __getitem__(self, i: Any) -> "IRanges":
        idx = positions_from_subscript(i, len(self), self._names)
        return self._take(idx)

    def _take(self, idx: IntArray) -> "IRanges":
        names = None if self._names is None else [self._names[k] for k in idx.tolist()]
        return IRanges._make(self._s()[idx], self._w()[idx], names, _mcols.take(self._mcols, idx))

    def concat(self, *others: "IRanges") -> "IRanges":
        """c(x, y, ...): names pad with "" when only some inputs are named; mcols union with NA fill."""
        parts = [self, *others]
        if any(not isinstance(p, IRanges) for p in parts):
            raise TypeError("can only concatenate IRanges with IRanges")
        start = np.concatenate([p._s() for p in parts])
        width = np.concatenate([p._w() for p in parts])
        names: Names | None = None
        if any(p._names is not None for p in parts):
            names = []
            for p in parts:
                names.extend(p._names if p._names is not None else [""] * len(p))
        mcols = _mcols.concat([(len(p), p._mcols) for p in parts])
        return IRanges._make(start, width, names, mcols)

    # ---- intra-range ------------------------------------------------------

    def shift(self, shift: Any) -> "IRanges":
        s = rep_len(as_int64(shift, "shift"), len(self), "shift")
        return IRanges._make(self._s() + s, self._w(), self._names, self._mcols)

    def narrow(self, start: Any = None, end: Any = None, width: Any = None) -> "IRanges":
        ns, nw = _intra.narrow_ranges(
            self._s(), self._w(), _int_arg(start, "start"), _int_arg(end, "end"), _int_arg(width, "width")
        )
        return IRanges._make(ns, nw, self._names, self._mcols)

    def resize(
        self, width: Any, fix: str | Sequence[str] | npt.NDArray[np.str_] = "start", use_names: bool = True
    ) -> "IRanges":
        fx = np.asarray([fix] if isinstance(fix, str) else list(fix))
        ns, nw = _intra.resize_ranges(self._s(), self._w(), as_int64(width, "width"), fx)
        return IRanges._make(ns, nw, self._names if use_names else None, self._mcols)

    def flank(
        self, width: Any, start: bool | Sequence[bool] | npt.NDArray[np.bool_] = True, both: bool = False
    ) -> "IRanges":
        at_start = as_bool_array(start if not isinstance(start, bool) else [start], len(self), "start")
        ns, nw = _intra.flank_ranges(self._s(), self._w(), as_int64(width, "width"), at_start, both)
        return IRanges._make(ns, nw, self._names, self._mcols)

    def promoters(self, upstream: Any = 2000, downstream: Any = 200) -> "IRanges":
        minus = np.zeros(len(self), dtype=bool)
        ns, nw = _intra.promoter_ranges(
            self._s(), self._w(), as_int64(upstream, "upstream"), as_int64(downstream, "downstream"), minus
        )
        return IRanges._make(ns, nw, self._names, self._mcols)

    # ---- inter-range --------------------------------------------------------

    def range(self, *others: "IRanges") -> "IRanges":
        both = self.concat(*others) if others else self
        s, w = _inter.range_of(both._s(), both._w())
        return IRanges._make(s, w)

    def reduce(
        self, drop_empty_ranges: bool = False, min_gapwidth: int = 1, with_revmap: bool = False
    ) -> "IRanges":
        s, w, rev = _inter.reduce_ranges(self._s(), self._w(), min_gapwidth, drop_empty_ranges, with_revmap)
        return IRanges._make(s, w, None, _revmap_frame(rev))

    def gaps(self, start: int | None = None, end: int | None = None) -> "IRanges":
        s, w = _inter.gaps_ranges(self._s(), self._w(), start, end)
        return IRanges._make(s, w)

    def disjoin(self, with_revmap: bool = False) -> "IRanges":
        s, w, rev = _inter.disjoin_ranges(self._s(), self._w(), with_revmap)
        return IRanges._make(s, w, None, _revmap_frame(rev))

    # ---- set operations -------------------------------------------------------

    def union(self, y: "IRanges") -> "IRanges":
        s, w = _inter.union_ranges(self._s(), self._w(), y._s(), y._w())
        return IRanges._make(s, w)

    def intersect(self, y: "IRanges") -> "IRanges":
        s, w = _inter.intersect_ranges(self._s(), self._w(), y._s(), y._w())
        return IRanges._make(s, w)

    def setdiff(self, y: "IRanges") -> "IRanges":
        s, w = _inter.setdiff_ranges(self._s(), self._w(), y._s(), y._w())
        return IRanges._make(s, w)

    def _check_parallel(self, y: "IRanges") -> None:
        if len(self) != len(y):
            raise ValueError("x and y must have the same length")

    def punion(self, y: "IRanges", fill_gap: bool = False) -> "IRanges":
        self._check_parallel(y)
        xs, xe, ys, ye = self._s(), self._e(), y._s(), y._e()
        if not fill_gap and np.any(np.maximum(xs, ys) > np.minimum(xe, ye) + 1):
            raise ValueError("some pair of ranges have a gap; use fill_gap=True to enforce their union")
        s = np.minimum(xs, ys)
        return IRanges._make(s, np.maximum(xe, ye) - s + 1, self._names)

    def pintersect(self, y: "IRanges", resolve_empty: str = "none") -> "IRanges":
        self._check_parallel(y)
        if resolve_empty not in ("none", "max.start", "start.x"):
            raise ValueError("resolve_empty must be 'none', 'max.start' or 'start.x'")
        xs, xe, ys, ye = self._s(), self._e(), y._s(), y._e()
        s = np.maximum(xs, ys)
        w = np.minimum(xe, ye) - s + 1
        ambiguous = w < 0
        if ambiguous.any():
            if resolve_empty == "none":
                raise ValueError("some intersections produce ambiguous empty ranges; use resolve_empty")
            if resolve_empty == "start.x":
                s = np.where(ambiguous, xs, s)
            w = np.where(ambiguous, 0, w)
        return IRanges._make(s, w, self._names)

    def psetdiff(self, y: "IRanges") -> "IRanges":
        self._check_parallel(y)
        xs, xe, ys, ye = self._s(), self._e(), y._s(), y._e()
        if np.any((ys > xs) & (ye < xe)):
            raise ValueError("some ranges in y have both end points strictly inside the range in x")
        hit = (ye - ys + 1 > 0) & (ys <= xe) & (ye >= xs)
        covers_all = hit & (ys <= xs) & (ye >= xe)
        left_cut = hit & (ys <= xs) & ~covers_all
        right_cut = hit & (ys > xs) & (ye >= xe)
        new_s = np.where(left_cut, ye + 1, xs)
        new_e = np.where(right_cut, ys - 1, xe)
        new_w = np.where(covers_all, 0, new_e - new_s + 1)
        new_s = np.where(covers_all, xs, new_s)
        return IRanges._make(new_s, new_w, self._names)

    # ---- overlaps ---------------------------------------------------------------

    def find_overlaps(
        self,
        subject: "IRanges",
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",
        select: str = "all",
    ) -> Hits | npt.NDArray[np.int64]:
        """Hits for select='all'; else an int64 array of subject indices with -1 for no hit."""
        if select not in ("all", "first", "last", "arbitrary"):
            raise ValueError("select must be 'all', 'first', 'last' or 'arbitrary'")
        q, s = _overlaps.overlap_pairs(
            self._s(), self._w(), subject._s(), subject._w(), maxgap, minoverlap, type
        )
        return _select_hits(q, s, len(self), len(subject), select)

    def count_overlaps(
        self, subject: "IRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any"
    ) -> npt.NDArray[np.int64]:
        hits = self.find_overlaps(subject, maxgap, minoverlap, type)
        assert isinstance(hits, Hits)
        return np.bincount(hits.query_hits, minlength=len(self)).astype(np.int64)

    def overlaps_any(
        self, subject: "IRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any"
    ) -> npt.NDArray[np.bool_]:
        return self.count_overlaps(subject, maxgap, minoverlap, type) > 0

    def subset_by_overlaps(
        self, ranges: "IRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any", invert: bool = False
    ) -> "IRanges":
        hit = self.overlaps_any(ranges, maxgap, minoverlap, type)
        return self[~hit if invert else hit]

    # ---- nearest ------------------------------------------------------------------

    def nearest(self, subject: "IRanges | None" = None, select: str = "arbitrary") -> Hits | npt.NDArray[np.int64]:
        """Closest subject (distance 0 if overlapping/adjacent). subject=None: nearest other range."""
        if select not in ("arbitrary", "all"):
            raise ValueError("select must be 'arbitrary' or 'all'")
        me = subject is None
        sj = self if subject is None else subject
        ties = _nearest.nearest_ties(self._s(), self._w(), sj._s(), sj._w(), None, me)
        return _ties_result(ties, len(self), len(sj), select)

    def precede(self, subject: "IRanges | None" = None, select: str = "first") -> Hits | npt.NDArray[np.int64]:
        if select not in ("first", "all"):
            raise ValueError("select must be 'first' or 'all'")
        sj = self if subject is None else subject  # R does not exclude self for precede/follow
        ties = _nearest.directional_ties(self._s(), self._w(), sj._s(), sj._w(), "precede")
        return _ties_result(ties, len(self), len(sj), select)

    def follow(self, subject: "IRanges | None" = None, select: str = "last") -> Hits | npt.NDArray[np.int64]:
        if select not in ("last", "all"):
            raise ValueError("select must be 'last' or 'all'")
        sj = self if subject is None else subject
        ties = _nearest.directional_ties(self._s(), self._w(), sj._s(), sj._w(), "follow")
        return _ties_result(ties, len(self), len(sj), select)

    def distance(self, other: "IRanges") -> npt.NDArray[np.int64]:
        return _nearest.distance_ranges(self._s(), self._w(), other._s(), other._w())

    # ---- coverage ---------------------------------------------------------------------

    def coverage(self, shift: Any = 0, width: int | None = None, weight: Any = 1) -> Rle:
        """Run-length coverage over positions 1..width (default: largest end, at least 0)."""
        sh = as_int64(shift, "shift")
        wt = np.asarray(weight)
        if wt.ndim == 0:
            wt = wt.reshape(1)
        if wt.dtype.kind not in "iuf":
            raise TypeError("weight must be numeric")
        if wt.dtype.kind == "f" and np.any(np.isnan(wt)):
            raise ValueError("weight cannot contain NAs")
        if len(wt) > max(len(self), 1):
            raise ValueError("weight is longer than x")
        if width is None:
            ends = self._e() + rep_len(sh, len(self), "shift") if len(self) else self._e()
            length = int(max(ends.max(), 0)) if len(self) else 0
        else:
            length = int(width)
        return _coverage.coverage_rle(self._s(), self._w(), length, sh, wt)

    # ---- ordering ---------------------------------------------------------------------

    def _keys(self) -> list[IntArray]:
        return [self._s(), self._w()]

    def order(self, decreasing: bool = False) -> IntArray:
        return _ordering.order_by(self._keys(), decreasing)

    def sort(self, decreasing: bool = False) -> "IRanges":
        return self._take(self.order(decreasing))

    def rank(self, ties_method: str = "average") -> npt.NDArray[np.generic]:
        return _ordering.rank_by(self._keys(), ties_method)

    def duplicated(self, from_last: bool = False) -> npt.NDArray[np.bool_]:
        return _ordering.duplicated_by(self._keys(), from_last)

    def unique(self, from_last: bool = False) -> "IRanges":
        return self[~self.duplicated(from_last)]

    # ---- misc -----------------------------------------------------------------------------

    def __repr__(self) -> str:
        shown = ", ".join(f"{s}-{s + w - 1}" for s, w in zip(self._start[:6].tolist(), self._width[:6].tolist()))
        more = ", ..." if len(self) > 6 else ""
        return f"IRanges(n={len(self)}: {shown}{more})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, IRanges):
            return NotImplemented
        same_cols = (self._mcols is None and other._mcols is None) or (
            self._mcols is not None and other._mcols is not None and self._mcols.equals(other._mcols)
        )
        return (
            np.array_equal(self._start, other._start)
            and np.array_equal(self._width, other._width)
            and self._names == other._names
            and same_cols
        )

    __hash__ = None  # type: ignore[assignment]


def _revmap_frame(revmap: list[IntArray] | None) -> pl.DataFrame | None:
    """revmap column: List(Int32) of 0-based input indices."""
    if revmap is None:
        return None
    series = pl.Series("revmap", [r.tolist() for r in revmap], dtype=pl.List(pl.Int32))
    return pl.DataFrame([series])


def _select_hits(
    q: IntArray, s: IntArray, nq: int, ns: int, select: str
) -> Hits | npt.NDArray[np.int64]:
    if select == "all":
        return Hits(q, s, nq, ns)
    out = np.full(nq, -1, dtype=np.int64)
    if len(q):
        if select == "last":
            np.maximum.at(out, q, s)
        else:  # first / arbitrary: lowest subject index
            big = np.iinfo(np.int64).max
            tmp = np.full(nq, big, dtype=np.int64)
            np.minimum.at(tmp, q, s)
            out = np.where(tmp == big, -1, tmp)
    return out


def _ties_result(
    ties: list[IntArray], nq: int, ns: int, select: str
) -> Hits | npt.NDArray[np.int64]:
    if select == "all":
        q = np.concatenate([np.full(len(t), i, dtype=np.int64) for i, t in enumerate(ties)] or [np.zeros(0, np.int64)])
        s = np.concatenate(ties or [np.zeros(0, np.int64)])
        return Hits(q, s, nq, ns)
    pick_last = select == "last"
    return np.array([(-1 if len(t) == 0 else int(t[-1] if pick_last else t[0])) for t in ties], dtype=np.int64)
