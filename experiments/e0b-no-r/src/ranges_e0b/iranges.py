"""IRanges: integer ranges, 1-based and closed (end = start + width - 1)."""

from collections.abc import Mapping, Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import coverage, inter, intra, nearest, ordering, overlaps, setops
from ._sew import (
    IntArray,
    concat_mcols,
    concat_names,
    empty_int,
    normalize_mcols,
    resolve_index,
    solve_sew,
    take_mcols,
)
from .overlaps import Hits
from .rle import Rle


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
        return self._mcols.clone()  # callers cannot mutate our metadata

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

    # --- intra-range ---------------------------------------------------------

    def shift(self, shift: Any) -> Self:
        """Move ranges by `shift` (scalar or per-range)."""
        return self._replace(*intra.shift(self._start, self._width, shift))

    def narrow(self, start: Any = None, end: Any = None, width: Any = None) -> Self:
        """Narrow relative to each range; negative start/end count from the range end."""
        return self._replace(*intra.narrow(self._start, self._width, start, end, width))

    def resize(self, width: Any, fix: str | Sequence[str] = "start") -> Self:
        """Set width, keeping `fix` ("start", "end" or "center") in place."""
        codes = intra.parse_fix(fix, len(self))  # type: ignore[arg-type]
        return self._replace(*intra.resize(self._start, self._width, width, codes))

    def flank(self, width: Any, start: Any = True, both: bool = False) -> Self:
        """Flanking ranges at the start (or end); negative width flanks inward."""
        at_start = intra.as_bool_array(start, len(self), "start")
        return self._replace(*intra.flank(self._start, self._width, width, at_start, both))

    def promoters(self, upstream: Any = 2000, downstream: Any = 200) -> Self:
        """`upstream` + `downstream` ranges around the start (all ranges treated as `+`)."""
        minus = np.zeros(len(self), dtype=bool)
        return self._replace(*intra.promoters(self._start, self._width, upstream, downstream, minus))

    # --- ordering ------------------------------------------------------------

    def _order_keys(self) -> list[IntArray]:
        return [self._start, self._width]

    def order(self, decreasing: bool = False) -> IntArray:
        """0-based stable permutation sorting by (start, width)."""
        return ordering.order_keys(self._order_keys(), decreasing)

    def sort(self, decreasing: bool = False) -> Self:
        return self._take(self.order(decreasing))

    def rank(self, ties_method: str = "average") -> Any:
        """1-based ranks by (start, width); float64 for "average", else int64."""
        return ordering.rank_keys(self._order_keys(), ties_method)

    def duplicated(self) -> Any:
        """Boolean mask; True for repeats of an earlier range (first occurrence False)."""
        return ordering.duplicated_keys(self._order_keys())

    def unique(self) -> Self:
        """Drop repeats, keeping first occurrences in original order."""
        return self[~self.duplicated()]

    # --- inter-range ---------------------------------------------------------
    # Results carry no names or mcols.

    def _group0(self) -> IntArray:
        return np.zeros(len(self), dtype=np.int64)

    @classmethod
    def _from_ends(cls, s: IntArray, e: IntArray) -> Self:
        return cls._new(s, e - s + 1)

    def range(self) -> Self:
        """Single range [min start, max end] (empty input -> empty)."""
        s, e, _ = inter.range_kernel(self._start, self.end, self._group0())
        return self._from_ends(s, e)

    def reduce(
        self,
        min_gapwidth: int = 1,
        drop_empty_ranges: bool = False,
        with_revmap: bool = False,
    ) -> Any:
        """Merge overlapping ranges and ranges separated by a gap < `min_gapwidth`.

        With `with_revmap=True` returns `(IRanges, list[ndarray])`; the lists hold
        0-based indices of the input ranges merged into each output range.
        """
        s, e, _, revmap = inter.reduce_kernel(
            self._start, self.end, self._group0(), min_gapwidth, drop_empty_ranges, with_revmap
        )
        out = self._from_ends(s, e)
        return (out, revmap) if with_revmap else out

    def gaps(self, start: int | None = None, end: int | None = None) -> Self:
        """Positions in [start, end] not covered; default bounds are min start / max end."""
        lo = start if start is not None else (int(self._start.min()) if len(self) else None)
        hi = end if end is not None else (int(self.end.max()) if len(self) else None)
        if lo is None or hi is None:
            return self._from_ends(empty_int(), empty_int())
        s, e, _ = inter.gaps_kernel(
            self._start,
            self.end,
            self._group0(),
            np.zeros(1, dtype=np.int64),
            np.array([lo], dtype=np.int64),
            np.array([hi], dtype=np.int64),
        )
        return self._from_ends(s, e)

    def disjoin(self) -> Self:
        """Disjoint pieces, split wherever any range starts or ends."""
        s, e, _ = inter.disjoin_kernel(self._start, self.end, self._group0())
        return self._from_ends(s, e)

    # --- set operations ------------------------------------------------------
    # Normalised results (reduced, non-empty), no names or mcols.

    def union(self, other: "IRanges") -> Self:
        s, e, _ = setops.union_kernel(
            self._start, self.end, self._group0(), other._start, other.end, other._group0()
        )
        return self._from_ends(s, e)

    def intersect(self, other: "IRanges") -> Self:
        s, e, _ = setops.sweep_kernel(
            self._start, self.end, self._group0(), other._start, other.end, other._group0(), "intersect"
        )
        return self._from_ends(s, e)

    def setdiff(self, other: "IRanges") -> Self:
        s, e, _ = setops.sweep_kernel(
            self._start, self.end, self._group0(), other._start, other.end, other._group0(), "setdiff"
        )
        return self._from_ends(s, e)

    # Parallel versions pair x[i] with y[i]; the result keeps x's names and mcols.

    def punion(self, other: "IRanges", fill_gap: bool = False) -> Self:
        s, e = setops.punion(self._start, self.end, other._start, other.end, fill_gap)
        return self._replace(s, e - s + 1)

    def pintersect(self, other: "IRanges", resolve_empty: str = "none") -> Self:
        s, e = setops.pintersect(self._start, self.end, other._start, other.end, resolve_empty)
        return self._replace(s, e - s + 1)

    def psetdiff(self, other: "IRanges") -> Self:
        s, e = setops.psetdiff(self._start, self.end, other._start, other.end)
        return self._replace(s, e - s + 1)

    # --- overlaps ------------------------------------------------------------

    def find_overlaps(
        self,
        subject: "IRanges | None" = None,
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",
        select: str = "all",
    ) -> Any:
        """Overlaps of each range (query) with `subject` (default: self).

        `select="all"` returns `Hits`; "first"/"last"/"arbitrary" return an int
        array with one 0-based subject index per query (`NO_HIT` = -1 if none).
        "arbitrary" is deterministic: same as "first".
        """
        subj = self if subject is None else subject
        qi, si = overlaps.overlap_pairs(self._start, self.end, subj._start, subj.end, type, maxgap, minoverlap)
        if select == "all":
            return Hits(qi, si, len(self), len(subj))
        return overlaps.select_hits(qi, si, len(self), select)

    def count_overlaps(self, subject: "IRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any") -> IntArray:
        hits = self.find_overlaps(subject, maxgap, minoverlap, type)
        return hits.count_per_query()

    def overlaps_any(self, subject: "IRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any") -> Any:
        return self.count_overlaps(subject, maxgap, minoverlap, type) > 0

    def subset_by_overlaps(
        self, subject: "IRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any", invert: bool = False
    ) -> Self:
        mask = self.overlaps_any(subject, maxgap, minoverlap, type)
        return self[~mask if invert else mask]

    # --- nearest -------------------------------------------------------------

    def _nearest_family(self, kind: str, subject: "IRanges | None", select: str, allowed: tuple[str, str]) -> Any:
        if select not in allowed:
            raise ValueError(f"select must be one of {allowed}, got {select!r}")
        subj = self if subject is None else subject
        zq, zs = np.zeros(len(self), np.int64), np.zeros(len(subj), np.int64)
        star_q, star_s = zq + nearest.STAR, zs + nearest.STAR
        qi, si = nearest.nearest_family(
            kind, self._start, self.end, zq, star_q, subj._start, subj.end, zs, star_s,
            ignore_strand=False, exclude_self=subject is None,
        )
        if select == "all":
            return Hits(qi, si, len(self), len(subj))
        return nearest.first_per_query(qi, si, len(self))

    def precede(self, subject: "IRanges | None" = None, select: str = "first") -> Any:
        """Subject with the smallest start > each range's end (adjacent counts).

        Returns 0-based indices (`NO_HIT` if none) or `Hits` for select="all".
        Ties go to the lowest index. `subject=None` means self, excluding self-matches.
        """
        return self._nearest_family("precede", subject, select, ("first", "all"))

    def follow(self, subject: "IRanges | None" = None, select: str = "first") -> Any:
        """Subject with the largest end < each range's start (adjacent counts)."""
        return self._nearest_family("follow", subject, select, ("first", "all"))

    def nearest(self, subject: "IRanges | None" = None, select: str = "arbitrary") -> Any:
        """Closest subject; overlapping or adjacent (gap 0) win; unchanged by the distance() redefinition.

        "arbitrary" is deterministic: the lowest index among ties.
        """
        return self._nearest_family("nearest", subject, select, ("arbitrary", "all"))

    def distance(self, other: "IRanges") -> IntArray:
        """Parallel distance; 0 if overlapping, else gap width + 1 (adjacent is 1)."""
        return nearest.distance(self._start, self.end, other._start, other.end)

    # --- coverage ------------------------------------------------------------

    def coverage(self, shift: int = 0, width: int | None = None, weight: Any = 1) -> Rle:
        """Per-position coverage over 1..width (default: max end after shift) as an `Rle`."""
        return coverage.coverage_kernel(self._start, self.end, shift, width, weight)
