"""IRanges: 1-based, closed integer ranges (``end = start + width - 1``)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _mcols, inter, intra, nearest, overlaps, setops
from .hits import Hits
from .rle import Rle
from . import coverage as _coverage
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

    # intra-range ---------------------------------------------------------
    def _with_coords(self, start: np.ndarray, width: np.ndarray, use_names: bool = True) -> Self:
        return self._make(start, width, self._names if use_names else None, self._mcols)

    def shift(self, shift: Any) -> Self:
        n = len(self)
        return self._with_coords(*intra.shift(self._start, self._width, intra.arg_vector(shift, n, "shift")))

    def narrow(self, start: Any = None, end: Any = None, width: Any = None) -> Self:
        n = len(self)
        a = [None if v is None else intra.arg_vector(v, n, nm) for v, nm in ((start, "start"), (end, "end"), (width, "width"))]
        return self._with_coords(*intra.narrow(self._start, self._width, *a))

    def resize(self, width: Any, fix: Any = "start", use_names: bool = True) -> Self:
        n = len(self)
        return self._with_coords(
            *intra.resize(self._start, self._width, intra.arg_vector(width, n, "width"), intra.str_vector(fix, n, "fix"), np.zeros(n, bool)),
            use_names,
        )

    def flank(self, width: Any, start: Any = True, both: bool = False, use_names: bool = True) -> Self:
        n = len(self)
        return self._with_coords(
            *intra.flank(
                self._start, self._width, intra.arg_vector(width, n, "width"),
                intra.flag_vector(start, n, "start"), both, np.zeros(n, bool),
            ),
            use_names,
        )

    def promoters(self, upstream: Any = 2000, downstream: Any = 200, use_names: bool = True) -> Self:
        n = len(self)
        return self._with_coords(
            *intra.promoters(
                self._start, self._width, intra.arg_vector(upstream, n, "upstream"),
                intra.arg_vector(downstream, n, "downstream"), np.zeros(n, bool),
            ),
            use_names,
        )

    # inter-range ---------------------------------------------------------
    @classmethod
    def _from_runs(
        cls, start: np.ndarray, width: np.ndarray, members: list[np.ndarray] | None = None
    ) -> Self:
        """Result of a collapsing op: no names; ``revmap`` mcols if members given."""
        mcols = None
        if members is not None:
            mcols = pl.DataFrame({"revmap": pl.Series("revmap", [m.tolist() for m in members], dtype=pl.List(pl.Int64))})
        return cls._make(start, width, None, mcols)

    def range(self, *others: Self) -> Self:
        allr = self.concat(self, *others) if others else self
        return self._from_runs(*inter.range_of(allr._start, allr._width))

    def reduce(self, drop_empty_ranges: bool = False, min_gapwidth: int = 1, with_revmap: bool = False) -> Self:
        s, w, members = inter.reduce_runs(self._start, self._width, int(min_gapwidth), drop_empty_ranges)
        return self._from_runs(s, w, members if with_revmap else None)

    def gaps(self, start: int | None = None, end: int | None = None) -> Self:
        for v, nm in ((start, "start"), (end, "end")):
            if v is not None and not isinstance(v, (int, np.integer)):
                raise TypeError(f"'{nm}' must be a single integer or None")
        return self._from_runs(*inter.gaps(self._start, self._width, start, end))

    def disjoin(self, with_revmap: bool = False) -> Self:
        s, w, members = inter.disjoin(self._start, self._width)
        return self._from_runs(s, w, members if with_revmap else None)

    # set operations -----------------------------------------------------
    def union(self, y: Self) -> Self:
        return self._from_runs(*setops.union(self._coords(), y._coords()))

    def intersect(self, y: Self) -> Self:
        return self._from_runs(*setops.intersect(self._coords(), y._coords()))

    def setdiff(self, y: Self) -> Self:
        return self._from_runs(*setops.setdiff(self._coords(), y._coords()))

    def _coords(self) -> setops.Coords:
        return self._start, self._width

    def _parallel(self, coords: setops.Coords) -> Self:
        """Parallel ops keep x's names and drop mcols (as R)."""
        return self._make(*coords, self._names)

    def punion(self, y: Self, fill_gap: bool = False) -> Self:
        return self._parallel(setops.punion(self._coords(), y._coords(), fill_gap))

    def pintersect(self, y: Self, resolve_empty: str = "none") -> Self:
        return self._parallel(setops.pintersect(self._coords(), y._coords(), resolve_empty))

    def psetdiff(self, y: Self) -> Self:
        return self._parallel(setops.psetdiff(self._coords(), y._coords()))

    # overlaps ------------------------------------------------------------
    def _overlap_pairs(self, subject: Self, maxgap: int, minoverlap: int, type: str) -> tuple[np.ndarray, np.ndarray]:
        return overlaps.overlap_pairs(
            self._start, self._width, subject._start, subject._width, int(maxgap), int(minoverlap), type
        )

    def find_overlaps(
        self, subject: Self | None = None, maxgap: int = -1, minoverlap: int = 0, type: str = "any", select: str = "all"
    ) -> Hits | np.ndarray:
        subject = self if subject is None else subject
        q, s = self._overlap_pairs(subject, maxgap, minoverlap, type)
        if select == "all":
            return Hits(q, s, len(self), len(subject))
        return overlaps.select_hits(q, s, len(self), select)

    def count_overlaps(self, subject: Self | None = None, maxgap: int = -1, minoverlap: int = 0, type: str = "any") -> np.ndarray:
        subject = self if subject is None else subject
        q, _ = self._overlap_pairs(subject, maxgap, minoverlap, type)
        return np.bincount(q, minlength=len(self)).astype(np.int64)

    def overlaps_any(self, subject: Self | None = None, maxgap: int = -1, minoverlap: int = 0, type: str = "any") -> np.ndarray:
        return self.count_overlaps(subject, maxgap, minoverlap, type) > 0

    def subset_by_overlaps(
        self, subject: Self | None = None, maxgap: int = -1, minoverlap: int = 0, type: str = "any", invert: bool = False
    ) -> Self:
        keep = self.overlaps_any(subject, maxgap, minoverlap, type)
        return self[~keep if invert else keep]

    # nearest -------------------------------------------------------------
    def _nearest_family(self, pairs_fn: Any, subject: Self, select: str, allowed: tuple[str, str]) -> Hits | np.ndarray:
        if select not in allowed:
            raise ValueError(f"'select' must be one of {allowed}")
        q, s = pairs_fn(self._start, self._width, subject._start, subject._width)
        if select == "all":
            return Hits(q, s, len(self), len(subject))
        return overlaps.select_hits(q, s, len(self), "last" if select == "last" else "first")

    def nearest(self, subject: Self, select: str = "arbitrary") -> Hits | np.ndarray:
        return self._nearest_family(nearest.nearest_pairs, subject, select, ("arbitrary", "all"))

    def precede(self, subject: Self, select: str = "first") -> Hits | np.ndarray:
        return self._nearest_family(nearest.precede_pairs, subject, select, ("first", "all"))

    def follow(self, subject: Self, select: str = "last") -> Hits | np.ndarray:
        return self._nearest_family(nearest.follow_pairs, subject, select, ("last", "all"))

    def distance(self, y: Self) -> np.ndarray:
        a, b = (self, y)
        if len(a) != len(b) and 1 not in (len(a), len(b)):
            raise ValueError("'x' and 'y' must have the same length (or one of length 1)")
        n = max(len(a), len(b)) if len(a) and len(b) else 0
        sa, wa = np.broadcast_to(a._start, (n,)), np.broadcast_to(a._width, (n,))
        sb, wb = np.broadcast_to(b._start, (n,)), np.broadcast_to(b._width, (n,))
        return nearest.distance(sa, wa, sb, wb)

    # coverage ------------------------------------------------------------
    def coverage(self, shift: Any = 0, width: int | None = None, weight: Any = 1) -> Rle:
        dense = _coverage.coverage_dense(
            self._start, self._width, shift, _coverage.check_length(width), weight
        )
        return Rle.encode(dense)

    def __repr__(self) -> str:
        rows = ", ".join(f"[{s},{e}]" for s, e in zip(self._start.tolist(), self.end.tolist(), strict=True))
        return f"IRanges({rows})"
