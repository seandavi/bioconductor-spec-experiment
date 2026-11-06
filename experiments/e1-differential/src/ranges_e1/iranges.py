"""IRanges: immutable vector of 1-based, closed integer ranges."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any, Self

import warnings

import numpy as np
import polars as pl

from . import _args, _meta
from ._args import flags
from ._kernels import coverage as coverage_kernel
from ._kernels import inter, intra, overlaps, setops
from ._kernels import order as order_kernel
from ._kernels import nearest as nearest_kernel
from .hits import Hits, overlap_result
from .rle import Rle
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
        if a.dtype.kind == "b":
            raise ValueError(f"'{what}' must be a numeric vector")
        if a.dtype.kind == "f" and not np.ma.isMaskedArray(v):
            if np.isnan(a).any():
                raise ValueError(f"'{what}' must be a numeric vector without NA's")
            a = np.trunc(a)
        elif a.dtype.kind not in "iufb" and not np.ma.isMaskedArray(v):
            raise ValueError(f"'{what}' must be a numeric vector")
        return _meta.recycle_arg(_meta.as_int_array(a if not np.ma.isMaskedArray(v) else v, f"'{what}'"), len(self), what, **kw)

    @flags("use_names")
    def shift(self, shift: Any = 0, use_names: bool = True) -> Self:
        by = self._int_arg(shift, "shift")
        s, w = intra.shift(self._start, self._width, by)
        return self._new(s, w, *self._meta_kept(use_names))

    @flags("use_names")
    def narrow(self, start: Any = None, end: Any = None, width: Any = None, use_names: bool = True) -> Self:
        """Narrow relative to each range; negative start/end count from the end (R)."""
        n = len(self)

        def sew(v: Any, what: str) -> np.ndarray:
            if v is None:
                return np.full(n, intra.NA, dtype=np.int64)
            if np.asarray(v).dtype.kind not in "iuf":
                raise ValueError(f"'{what}' must be a vector of integers")
            if np.ma.isMaskedArray(v):
                arr = np.where(np.ma.getmaskarray(v), intra.NA, np.ma.getdata(v)).astype(np.int64)
            else:
                f = np.atleast_1d(np.asarray(v, dtype=float))
                arr = np.where(np.isnan(f), intra.NA, np.trunc(np.nan_to_num(f))).astype(np.int64)
            return _meta.recycle_arg(arr, n, what, warn=False)

        s, w = intra.narrow(self._start, self._width, sew(start, "start"), sew(end, "end"), sew(width, "width"))
        return self._new(s, w, *self._meta_kept(use_names))

    @flags("use_names")
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

    @flags("both", "use_names")
    def flank(self, width: Any, start: Any = True, both: bool = False, use_names: bool = True) -> Self:
        by = self._int_arg(width, "width")
        flags = np.asarray(start)
        if flags.dtype.kind != "b":
            raise ValueError("'start' must be logical without NA's")
        at_start = _meta.recycle_arg(np.atleast_1d(flags), len(self), "start", longer="truncate")
        s, w = intra.flank(self._start, self._width, by, at_start, bool(both))
        return self._new(s, w, *self._meta_kept(use_names))

    @flags("use_names")
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

    @flags("drop_empty_ranges", "with_revmap")
    def reduce(
        self, drop_empty_ranges: bool = False, min_gapwidth: int = 1, with_revmap: bool = False
    ) -> Self:
        gap = _args.scalar_int(min_gapwidth, "min.gapwidth")
        assert gap is not None
        if gap < 0:
            raise ValueError("'min.gapwidth' must be non-negative")
        red = inter.reduce(self._group0(), self._start, self._width, gap, drop_empty_ranges)
        mc = _meta.revmap_frame(red.members) if with_revmap else None
        return self._new(red.start, red.width, None, mc)

    def gaps(self, start: int | None = None, end: int | None = None) -> Self:
        """Uncovered stretches within [start, end]; defaults to the extent of non-empty ranges."""
        start, end = _args.scalar_int(start, "start", allow_na=True), _args.scalar_int(end, "end", allow_na=True)
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

    @flags("with_revmap")
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

    @flags("fill_gap")
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

    @flags("drop_self", "drop_redundant")
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
        type, select, gap, minov = overlaps.normalize_args(type, select, maxgap, minoverlap)
        subj = self if subject is None else subject
        q, s = self._pairs(subj, gap, minov, type)
        return overlap_result(
            q, s, len(self), len(subj), select=select, self_hits=subject is None,
            drop_self=drop_self, drop_redundant=drop_redundant,
            query_is_shorter=overlaps.uses_query_nclist(len(self), len(subj)),
        )

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

    @flags("invert")
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

    # -- nearest / precede / follow / distance -------------------------------
    def _pick(self, q: np.ndarray, s: np.ndarray) -> np.ma.MaskedArray:
        """Per-query subject index from (q, s) pairs with at most one pair per query."""
        out = np.zeros(len(self), dtype=np.int64)
        has = np.zeros(len(self), dtype=bool)
        out[q], has[q] = s, True
        return np.ma.masked_array(out, mask=~has)

    def precede(self, subject: IRanges | None = None, select: str = "first") -> Hits | np.ma.MaskedArray:
        """Subject(s) starting right after each range ends (smallest start > end)."""
        select = _args.match_arg(select, ("first", "all"))
        subj = self if subject is None else subject
        q, s = nearest_kernel.precede(self.end, subj._start, select == "all")
        return Hits(q, s, len(self), len(subj)) if select == "all" else self._pick(q, s)

    def follow(self, subject: IRanges | None = None, select: str = "last") -> Hits | np.ma.MaskedArray:
        """Subject(s) ending right before each range starts (largest end < start)."""
        select = _args.match_arg(select, ("last", "all"))
        subj = self if subject is None else subject
        q, s = nearest_kernel.follow(self._start, subj.end, select == "all")
        return Hits(q, s, len(self), len(subj)) if select == "all" else self._pick(q, s)

    def nearest(self, subject: IRanges | None = None, select: str = "arbitrary") -> Hits | np.ma.MaskedArray:
        """Nearest subject range: overlapping/adjacent first, else the closer of precede/follow.

        Ties between the two sides go to the follow side (select="arbitrary") or keep both
        (select="all"), as R does.
        """
        select = _args.match_arg(select, ("arbitrary", "all"))
        subj = self if subject is None else subject
        n = len(self)
        ol = self.find_overlaps(subject, maxgap=0, select=select, drop_self=subject is None)
        olv = ol.select("first") if isinstance(ol, Hits) else ol
        miss = np.flatnonzero(np.ma.getmaskarray(olv))
        sub = self._take(miss)
        all_ = select == "all"
        bq, bs = nearest_kernel.precede(sub.end, subj._start, all_)
        aq, as_ = nearest_kernel.follow(sub._start, subj.end, all_)
        # first subject on each side, per missing query (NA -> -1 here, internal only)
        def first_of(q: np.ndarray, s: np.ndarray) -> np.ndarray:
            out = np.full(len(miss), -1, dtype=np.int64)
            out[q[::-1]] = s[::-1]
            return out
        before, after = first_of(bq, bs), first_of(aq, as_)
        has_b, has_a = before >= 0, after >= 0
        big = np.iinfo(np.int64).max // 4
        if len(subj) == 0:
            before_dist = after_dist = np.full(len(miss), big, dtype=np.int64)
        else:
            before_dist = np.where(has_b, subj._start[np.maximum(before, 0)] - sub.end, big)
            after_dist = np.where(has_a, sub._start - subj.end[np.maximum(after, 0)], big)
        both = has_b & has_a
        left = np.where(both, before_dist < after_dist, ~has_a)
        if not all_:
            res = np.ma.getdata(olv).copy()
            res[miss] = np.where(left, before, after)
            mask = np.ma.getmaskarray(olv).copy()
            mask[miss] = ~np.where(left, has_b, has_a)
            return np.ma.masked_array(res, mask=mask)
        right = ~left
        tie = both & (before_dist == after_dist)
        left = left | tie
        assert isinstance(ol, Hits)
        keep_b, keep_a = left[bq], right[aq]
        qs_ = np.concatenate([ol.query_hits, miss[bq[keep_b]], miss[aq[keep_a]]])
        ss_ = np.concatenate([ol.subject_hits, bs[keep_b], as_[keep_a]])
        order = np.lexsort((ss_, qs_))
        return Hits(qs_[order], ss_[order], n, len(subj))

    def distance(self, other: IRanges) -> np.ma.MaskedArray:
        """Parallel gap between ranges: 0 when they overlap or touch (R: distance).

        A gap beyond int32 is NA (masked), as R's integer arithmetic overflows to NA.
        """
        a, b = len(self), len(other)
        if a == 0 or b == 0:
            return np.ma.masked_array(np.empty(0, dtype=np.int64))
        n = max(a, b)
        if n % a or n % b:
            warnings.warn("longer object length is not a multiple of shorter object length", UserWarning, stacklevel=2)
        xs, xe = np.resize(self._start, n), np.resize(self.end, n)
        ys, ye = np.resize(other._start, n), np.resize(other.end, n)
        d = np.maximum(np.maximum(xs, ys) - np.minimum(xe, ye) - 1, 0)
        over = d > INT_MAX
        if over.any():
            warnings.warn("NAs produced by integer overflow", UserWarning, stacklevel=2)
        return np.ma.masked_array(np.where(over, 0, d), mask=over)

    # -- coverage ------------------------------------------------------------
    def _coverage_inputs(
        self, shift: Any, weight: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        """Shifted starts and weights, recycled to len(self) (warn if not a divisor)."""
        n = len(self)
        sh = self._int_arg(shift, "shift", longer="truncate")
        if np.ma.isMaskedArray(weight):
            raise ValueError("'weight' cannot contain NAs")
        wt = np.atleast_1d(np.asarray(weight))
        if wt.dtype.kind not in "iuf":
            raise ValueError("'weight' must be an integer or numeric vector")
        if wt.dtype.kind == "f" and np.isnan(wt).any():
            raise ValueError("'weight' cannot contain NAs")
        wt = wt.astype(np.int64) if wt.dtype.kind in "iu" else wt.astype(np.float64)
        wt = _meta.recycle_arg(wt, n, "weight", longer="truncate")
        return self._start + sh, wt

    def coverage(self, shift: Any = 0, width: int | None = None, weight: Any = 1) -> Rle:
        """Depth of coverage over positions 1..width as an Rle.

        Default width is the largest end after shifting (0 if none).
        """
        start, wt = self._coverage_inputs(shift, weight)
        w = _args.scalar_int(width, "width", allow_na=True) if width is not None else None
        if w is None:
            total = int(max((start + self._width - 1).max(initial=0), 0))
        else:
            total = w
            if total < 0:
                raise ValueError("'width' cannot be negative")
        vals, lens = coverage_kernel.coverage(start, self._width, wt, total)
        return Rle(vals, lens)

    # -- ordering ------------------------------------------------------------
    def _order_keys(self) -> list[np.ndarray]:
        """Sort keys, most significant first: start, then width."""
        return [self._start, self._width]

    @staticmethod
    def _flag(v: Any, what: str) -> bool:
        if not isinstance(v, (bool, np.bool_)):
            raise ValueError(f"'{what}' must be TRUE or FALSE")
        return bool(v)

    @flags("decreasing")
    def order(self, decreasing: bool = False) -> np.ndarray:
        """0-based permutation sorting by (start, width); ties keep input order."""
        return order_kernel.order(self._order_keys(), decreasing)

    @flags("decreasing")
    def sort(self, decreasing: bool = False) -> Self:
        return self._take(self.order(decreasing))

    def rank(self, ties_method: str = "average") -> np.ndarray:
        """1-based ranks by (start, width); float for "average", int otherwise.

        R's own default fails to match its ties.method choices; Python defaults to "average" (ADR-0004).
        """
        return order_kernel.rank(self._order_keys(), _args.match_arg(ties_method, order_kernel.TIES, "ties.method"))

    @flags("from_last")
    def duplicated(self, from_last: bool = False) -> np.ndarray:
        return order_kernel.duplicated(self._order_keys(), from_last)

    @flags("from_last")
    def unique(self, from_last: bool = False) -> Self:
        return self[~self.duplicated(from_last)]
