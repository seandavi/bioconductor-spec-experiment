"""GRanges: IRanges + seqnames + strand + seqinfo."""

import warnings
from collections.abc import Mapping, Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import coverage, inter, intra, nearest, ordering, overlaps, setops
from ._sew import (
    NO_HIT,
    BoolArray,
    IntArray,
    concat_mcols,
    empty_int,
    normalize_mcols,
    resolve_index,
    take_mcols,
)
from .iranges import IRanges
from .overlaps import Hits
from .rle import Rle
from .seqinfo import Seqinfo

STRANDS = ("+", "-", "*")
PLUS, MINUS, STAR = 0, 1, 2
_STRAND_CODE = {"+": PLUS, "-": MINUS, "*": STAR}


def _remap(codes: IntArray, old: Seqinfo, new: Seqinfo) -> IntArray:
    if old is new or old.seqnames == new.seqnames:
        return codes
    mapping = np.array([new.index(n) for n in old.seqnames], dtype=np.int64)
    return mapping[codes] if len(codes) else codes


def _strand_codes(strand: Any, n: int) -> IntArray:
    vals = [strand] if isinstance(strand, str) else list(strand)
    try:
        codes = np.array([_STRAND_CODE[v] for v in vals], dtype=np.int64)
    except KeyError as err:
        raise ValueError(f"strand must be '+', '-' or '*', got {err.args[0]!r}") from None
    if len(codes) == 1:
        return np.repeat(codes, n)
    if len(codes) != n:
        raise ValueError(f"strand has length {len(codes)}, expected {n}")
    return codes


class GRanges:
    """Immutable genomic ranges: IRanges plus per-range seqname and strand and a Seqinfo.

    `seqnames` are held as codes into `seqinfo.seqnames`, so the seqlevels order
    (not alphabetical order) drives sorting, as in R.
    """

    __slots__ = ("_seq", "_strand", "_ranges", "_mcols", "_seqinfo")

    def __init__(
        self,
        seqnames: str | Sequence[str],
        ranges: IRanges,
        strand: str | Sequence[str] = "*",
        mcols: pl.DataFrame | Mapping[str, Any] | None = None,
        seqinfo: Seqinfo | None = None,
    ) -> None:
        n = len(ranges)
        names = [seqnames] * n if isinstance(seqnames, str) else list(seqnames)
        if len(names) != n:
            raise ValueError(f"seqnames has length {len(names)}, expected {n}")
        if seqinfo is None:
            seqinfo = Seqinfo(list(dict.fromkeys(names)))
        index = {nm: i for i, nm in enumerate(seqinfo.seqnames)}
        try:
            seq = np.array([index[nm] for nm in names], dtype=np.int64)
        except KeyError as err:
            raise ValueError(f"seqname {err.args[0]!r} not in seqinfo") from None
        self._init(seq, _strand_codes(strand, n), ranges.with_mcols(None), normalize_mcols(mcols, n), seqinfo)

    def _init(
        self, seq: IntArray, strand: IntArray, ranges: IRanges, mcols: pl.DataFrame, seqinfo: Seqinfo, check: bool = True
    ) -> None:
        seq.setflags(write=False)
        strand.setflags(write=False)
        self._seq, self._strand, self._ranges, self._mcols, self._seqinfo = seq, strand, ranges, mcols, seqinfo
        if check:
            self._warn_out_of_bound()

    @classmethod
    def _make(
        cls,
        seq: IntArray,
        strand: IntArray,
        ranges: IRanges,
        seqinfo: Seqinfo,
        mcols: pl.DataFrame | None = None,
        check: bool = True,
    ) -> Self:
        obj = object.__new__(cls)
        obj._init(np.ascontiguousarray(seq, dtype=np.int64), np.ascontiguousarray(strand, dtype=np.int64),
                  ranges, pl.DataFrame() if mcols is None else mcols, seqinfo, check)
        return obj

    # --- bounds ----------------------------------------------------------------

    def _bounds(self) -> tuple[BoolArray, IntArray]:
        """(applicable, seqlength) per range: known length and not circular."""
        lengths = np.array([-1 if v is None else v for v in self._seqinfo.seqlengths] or [-1], dtype=np.int64)
        circ = np.array([bool(v) for v in self._seqinfo.is_circular] or [False])
        if len(self._seq) == 0:
            return np.zeros(0, dtype=bool), empty_int()
        ln = lengths[self._seq]
        return (ln >= 0) & ~circ[self._seq], ln

    def _out_of_bound(self) -> BoolArray:
        applicable, ln = self._bounds()
        return applicable & ((self._ranges.start < 1) | (self._ranges.end > ln))

    def _warn_out_of_bound(self) -> None:
        bad = self._out_of_bound()
        if bad.any():
            levels = sorted({self._seqinfo.seqnames[c] for c in self._seq[bad]})
            warnings.warn(
                f"GRanges contains {int(bad.sum())} out-of-bound range(s) on sequence(s) "
                f"{', '.join(levels)}; use trim() to trim them",
                stacklevel=3,
            )

    def trim(self) -> Self:
        """Clip ranges to [1, seqlength]; fully outside ranges become empty at the boundary."""
        applicable, ln = self._bounds()
        s, e = self._ranges.start, self._ranges.end
        s2 = np.where(applicable, np.maximum(s, 1), s)
        e2 = np.where(applicable, np.minimum(e, ln), e)
        s2 = np.where(applicable, np.minimum(s2, ln + 1), s2)
        e2 = np.where(applicable, np.maximum(e2, 0), e2)
        e2 = np.maximum(e2, s2 - 1)
        return self._with_coords(s2, e2 - s2 + 1)

    # --- accessors -------------------------------------------------------------

    @property
    def seqnames(self) -> list[str]:
        names = np.array(self._seqinfo.seqnames, dtype=object)
        return names[self._seq].tolist() if len(self._seq) else []

    @property
    def strand(self) -> list[str]:
        return [STRANDS[c] for c in self._strand]

    @property
    def ranges(self) -> IRanges:
        return self._ranges

    @property
    def start(self) -> IntArray:
        return self._ranges.start

    @property
    def end(self) -> IntArray:
        return self._ranges.end

    @property
    def width(self) -> IntArray:
        return self._ranges.width

    @property
    def names(self) -> list[str] | None:
        return self._ranges.names

    @property
    def mcols(self) -> pl.DataFrame:
        return self._mcols.clone()  # callers cannot mutate our metadata

    @property
    def seqinfo(self) -> Seqinfo:
        return self._seqinfo

    @property
    def seqlevels(self) -> list[str]:
        return list(self._seqinfo.seqnames)

    @property
    def seqlengths(self) -> dict[str, int | None]:
        return dict(zip(self._seqinfo.seqnames, self._seqinfo.seqlengths))

    @property
    def is_circular(self) -> dict[str, bool | None]:
        return dict(zip(self._seqinfo.seqnames, self._seqinfo.is_circular))

    @property
    def genome(self) -> dict[str, str | None]:
        return dict(zip(self._seqinfo.seqnames, self._seqinfo.genome))

    def __len__(self) -> int:
        return len(self._seq)

    def __repr__(self) -> str:
        n = len(self)
        head = f"GRanges of length {n}, {len(self._seqinfo)} seqlevels"
        if self._mcols.width:
            head += f", mcols ({', '.join(self._mcols.columns)})"
        rows = [head]
        shown = list(range(n)) if n <= 10 else list(range(5)) + [-1] + list(range(n - 5, n))
        sn = self.seqnames if n <= 10 else None
        for i in shown:
            if i == -1:
                rows.append("  ...")
                continue
            name = sn[i] if sn is not None else self._seqinfo.seqnames[self._seq[i]]
            rows.append(f"  [{i}] {name}:{self.start[i]}-{self.end[i]} [{STRANDS[self._strand[i]]}]")
        return "\n".join(rows)

    def with_names(self, names: Sequence[str] | None) -> Self:
        return self._make(self._seq, self._strand, self._ranges.with_names(names), self._seqinfo, self._mcols, check=False)

    def with_mcols(self, mcols: pl.DataFrame | Mapping[str, Any] | None) -> Self:
        return self._make(self._seq, self._strand, self._ranges, self._seqinfo, normalize_mcols(mcols, len(self)), check=False)

    def with_seqinfo(self, seqinfo: Seqinfo) -> Self:
        """Replace seqinfo; every used seqname must exist in the new one (codes are remapped)."""
        used = {self._seqinfo.seqnames[c] for c in np.unique(self._seq)}
        missing = used - set(seqinfo.seqnames)
        if missing:
            raise ValueError(f"new seqinfo lacks used seqlevels: {sorted(missing)}")
        return self._make(_remap(self._seq, self._seqinfo, seqinfo), self._strand, self._ranges, seqinfo, self._mcols)

    def equals(self, other: "GRanges") -> bool:
        return (
            len(self) == len(other)
            and self._seqinfo == other._seqinfo
            and bool(np.array_equal(self._seq, other._seq))
            and bool(np.array_equal(self._strand, other._strand))
            and self._ranges.equals(other._ranges)
            and self._mcols.equals(other._mcols)
        )

    def _with_coords(self, start: IntArray, width: IntArray) -> Self:
        """Same seq/strand/names/mcols/seqinfo, new coordinates."""
        return self._make(self._seq, self._strand, self._ranges._replace(start, width), self._seqinfo, self._mcols)

    # --- subsetting and combining ------------------------------------------------

    def __getitem__(self, i: Any) -> Self:
        return self._take(resolve_index(i, len(self), self._ranges._names))

    def _take(self, idx: IntArray) -> Self:
        return self._make(
            self._seq[idx], self._strand[idx], self._ranges._take(idx), self._seqinfo,
            take_mcols(self._mcols, idx), check=False,
        )

    @classmethod
    def combine(cls, items: Sequence["GRanges"]) -> Self:
        """Concatenate; seqinfo is merged (conflicts raise); mcols columns must match."""
        if not items:
            return cls("*", IRanges())  # type: ignore[arg-type]
        info = items[0]._seqinfo
        for x in items[1:]:
            if x._seqinfo != info:
                info = info.merge(x._seqinfo)
        seq = np.concatenate([_remap(x._seq, x._seqinfo, info) for x in items])
        ranges = IRanges.combine([x._ranges for x in items])
        mcols = concat_mcols([x._mcols for x in items])
        return cls._make(seq, np.concatenate([x._strand for x in items]), ranges, info, mcols)

    def concat(self, *others: "GRanges") -> Self:
        return self.combine([self, *others])

    # --- intra-range -----------------------------------------------------------
    # shift and narrow ignore strand; resize, flank and promoters treat `-` as
    # reversed (their "start" is the range end) and `*` as `+`.

    def _minus(self, ignore_strand: bool = False) -> BoolArray:
        if ignore_strand:
            return np.zeros(len(self), dtype=bool)
        return self._strand == MINUS

    def shift(self, shift: Any) -> Self:
        return self._with_coords(*intra.shift(self._ranges.start, self._ranges.width, shift))

    def narrow(self, start: Any = None, end: Any = None, width: Any = None) -> Self:
        return self._with_coords(*intra.narrow(self._ranges.start, self._ranges.width, start, end, width))

    def resize(self, width: Any, fix: str | Sequence[str] = "start", ignore_strand: bool = False) -> Self:
        codes = intra.parse_fix(fix, len(self))  # type: ignore[arg-type]
        swap = self._minus(ignore_strand) & (codes != intra.FIX_CENTER)
        codes = np.where(swap, 1 - codes, codes)
        return self._with_coords(*intra.resize(self._ranges.start, self._ranges.width, width, codes))

    def flank(self, width: Any, start: Any = True, both: bool = False, ignore_strand: bool = False) -> Self:
        at_start = intra.as_bool_array(start, len(self), "start") != self._minus(ignore_strand)
        return self._with_coords(*intra.flank(self._ranges.start, self._ranges.width, width, at_start, both))

    def promoters(self, upstream: Any = 2000, downstream: Any = 200, ignore_strand: bool = False) -> Self:
        return self._with_coords(
            *intra.promoters(
                self._ranges.start, self._ranges.width, upstream, downstream, self._minus(ignore_strand)
            )
        )

    # --- ordering --------------------------------------------------------------
    # Key: (seqlevel position, strand + < - < *, start, width).

    def _order_keys(self) -> list[IntArray]:
        return [self._seq, self._strand, self._ranges.start, self._ranges.width]

    def order(self, decreasing: bool = False) -> IntArray:
        return ordering.order_keys(self._order_keys(), decreasing)

    def sort(self, decreasing: bool = False) -> Self:
        return self._take(self.order(decreasing))

    def rank(self, ties_method: str = "average") -> Any:
        return ordering.rank_keys(self._order_keys(), ties_method)

    def duplicated(self) -> Any:
        return ordering.duplicated_keys(self._order_keys())

    def unique(self) -> Self:
        return self[~self.duplicated()]

    # --- inter-range -------------------------------------------------------------
    # Grouped by (seqlevel, strand); results are sorted by seqlevel, strand
    # (+ - *), start, and carry no names or mcols. ignore_strand=True puts
    # every range in the `*` group.

    def _groups(self, ignore_strand: bool) -> IntArray:
        return self._seq * 3 + (STAR if ignore_strand else self._strand)

    def _from_groups(self, s: IntArray, e: IntArray, g: IntArray, seqinfo: Seqinfo | None = None) -> Self:
        return self._make(g // 3, g % 3, IRanges._from_ends(s, e), seqinfo or self._seqinfo, check=False)

    def range(self, ignore_strand: bool = False) -> Self:
        s, e, g = inter.range_kernel(self._ranges.start, self._ranges.end, self._groups(ignore_strand))
        return self._from_groups(s, e, g)

    def reduce(
        self,
        min_gapwidth: int = 1,
        drop_empty_ranges: bool = False,
        with_revmap: bool = False,
        ignore_strand: bool = False,
    ) -> Any:
        """As IRanges.reduce, per (seqlevel, strand); with_revmap returns `(GRanges, list)`."""
        s, e, g, revmap = inter.reduce_kernel(
            self._ranges.start, self._ranges.end, self._groups(ignore_strand), min_gapwidth, drop_empty_ranges, with_revmap
        )
        out = self._from_groups(s, e, g)
        return (out, revmap) if with_revmap else out

    def disjoin(self, ignore_strand: bool = False) -> Self:
        s, e, g = inter.disjoin_kernel(self._ranges.start, self._ranges.end, self._groups(ignore_strand))
        return self._from_groups(s, e, g)

    def gaps(
        self,
        start: int | Mapping[str, int] = 1,
        end: int | Mapping[str, int] | None = None,
        ignore_strand: bool = False,
    ) -> Self:
        """Uncovered positions within [start, end] for every seqlevel x strand.

        Every seqlevel yields all three strands (only `*` with ignore_strand), so
        a strand without ranges gives the full-length gap. `end` defaults to the
        seqlength, else the largest end on that seqlevel; a seqlevel with neither
        is skipped.
        """
        names = self._seqinfo.seqnames
        max_end: dict[int, int] = {}
        if len(self):
            for k in np.unique(self._seq):
                max_end[int(k)] = int(self._ranges.end[self._seq == k].max())
        strands = (STAR,) if ignore_strand else (PLUS, MINUS, STAR)
        groups: list[int] = []
        los: list[int] = []
        his: list[int] = []
        for k, name in enumerate(names):
            lo = start[name] if isinstance(start, Mapping) else start
            hi: int | None
            if end is None:
                hi = self._seqinfo.seqlengths[k]
                if hi is None:
                    hi = max_end.get(k)
            else:
                hi = end[name] if isinstance(end, Mapping) else end
            if hi is None:
                continue
            for st in strands:
                groups.append(k * 3 + st)
                los.append(lo)
                his.append(hi)
        s, e, g = inter.gaps_kernel(
            self._ranges.start, self._ranges.end, self._groups(ignore_strand),
            np.array(groups, dtype=np.int64), np.array(los, dtype=np.int64), np.array(his, dtype=np.int64),
        )
        return self._from_groups(s, e, g)

    # --- set operations ----------------------------------------------------------

    def _align(self, other: "GRanges") -> tuple[Seqinfo, IntArray, IntArray]:
        """Merged seqinfo and both objects' seq codes expressed in it."""
        info = self._seqinfo if other._seqinfo == self._seqinfo else self._seqinfo.merge(other._seqinfo)
        return info, _remap(self._seq, self._seqinfo, info), _remap(other._seq, other._seqinfo, info)

    def _set_args(
        self, other: "GRanges", ignore_strand: bool
    ) -> tuple[Seqinfo, tuple[IntArray, IntArray, IntArray, IntArray, IntArray, IntArray]]:
        info, a_seq, b_seq = self._align(other)
        star = STAR if ignore_strand else None
        ga = a_seq * 3 + (star if star is not None else self._strand)
        gb = b_seq * 3 + (star if star is not None else other._strand)
        return info, (self._ranges.start, self._ranges.end, ga, other._ranges.start, other._ranges.end, gb)

    def union(self, other: "GRanges", ignore_strand: bool = False) -> Self:
        info, args = self._set_args(other, ignore_strand)
        s, e, g = setops.union_kernel(*args)
        return self._from_groups(s, e, g, info)

    def intersect(self, other: "GRanges", ignore_strand: bool = False) -> Self:
        info, args = self._set_args(other, ignore_strand)
        s, e, g = setops.sweep_kernel(*args, "intersect")
        return self._from_groups(s, e, g, info)

    def setdiff(self, other: "GRanges", ignore_strand: bool = False) -> Self:
        info, args = self._set_args(other, ignore_strand)
        s, e, g = setops.sweep_kernel(*args, "setdiff")
        return self._from_groups(s, e, g, info)

    def _parallel_args(
        self, other: "GRanges", ignore_strand: bool, keep_x_strand: bool = False
    ) -> tuple[Seqinfo, IntArray, IntArray]:
        """Merged seqinfo, seq codes and result strands for pairwise ops."""
        if len(self) != len(other):
            raise ValueError("parallel operations require equal lengths")
        info, a_seq, b_seq = self._align(other)
        if np.any(a_seq != b_seq):
            raise ValueError("parallel operations require identical seqnames")
        a, b = self._strand, other._strand
        if ignore_strand:
            return info, a_seq, a.copy() if keep_x_strand else np.full(len(a), STAR, dtype=np.int64)
        if np.any((a != b) & (a != STAR) & (b != STAR)):
            raise ValueError("parallel operations require compatible strands (equal, or `*`)")
        return info, a_seq, a if keep_x_strand else np.where(a == STAR, b, a)

    def _parallel_result(self, info: Seqinfo, seq: IntArray, strand: IntArray, s: IntArray, e: IntArray) -> Self:
        return self._make(seq, strand, self._ranges._replace(s, e - s + 1), info, self._mcols)

    def punion(self, other: "GRanges", fill_gap: bool = False, ignore_strand: bool = False) -> Self:
        info, seq, strand = self._parallel_args(other, ignore_strand)
        s, e = setops.punion(self._ranges.start, self._ranges.end, other._ranges.start, other._ranges.end, fill_gap)
        return self._parallel_result(info, seq, strand, s, e)

    def pintersect(self, other: "GRanges", resolve_empty: str = "none", ignore_strand: bool = False) -> Self:
        info, seq, strand = self._parallel_args(other, ignore_strand)
        s, e = setops.pintersect(
            self._ranges.start, self._ranges.end, other._ranges.start, other._ranges.end, resolve_empty
        )
        return self._parallel_result(info, seq, strand, s, e)

    def psetdiff(self, other: "GRanges", ignore_strand: bool = False) -> Self:
        info, seq, strand = self._parallel_args(other, ignore_strand, keep_x_strand=True)
        s, e = setops.psetdiff(self._ranges.start, self._ranges.end, other._ranges.start, other._ranges.end)
        return self._parallel_result(info, seq, strand, s, e)

    # --- overlaps ----------------------------------------------------------------
    # Pairs must share a seqlevel; strands must be equal or either `*` unless
    # ignore_strand. For type "start"/"end", a pair involving a `-` range compares
    # the opposite ends (5'/3' semantics). Circular wraparound is not handled.

    def _overlap_pairs(
        self, subj: "GRanges", maxgap: int, minoverlap: int, type: str, ignore_strand: bool
    ) -> tuple[IntArray, IntArray]:
        overlaps.check_args(type, maxgap, minoverlap)
        _, q_seq, s_seq = self._align(subj)
        qs, qe, ss, se = self._ranges.start, self._ranges.end, subj._ranges.start, subj._ranges.end
        flips = type in ("start", "end") and not ignore_strand
        other_type = {"start": "end", "end": "start"}.get(type, type)
        out_q: list[IntArray] = []
        out_s: list[IntArray] = []
        for k in np.unique(q_seq):
            qidx, sidx = np.flatnonzero(q_seq == k), np.flatnonzero(s_seq == k)
            if len(qidx) == 0 or len(sidx) == 0:
                continue
            args = (qs[qidx], qe[qidx], ss[sidx], se[sidx])
            lq, ls = overlaps.overlap_pairs(*args, type, maxgap, minoverlap)
            q, s = qidx[lq], sidx[ls]
            qstr, sstr = self._strand[q], subj._strand[s]
            if not ignore_strand:
                ok = (qstr == STAR) | (sstr == STAR) | (qstr == sstr)
                q, s, qstr, sstr = q[ok], s[ok], qstr[ok], sstr[ok]
            if flips:
                minus = (qstr == MINUS) | (sstr == MINUS)
                lq2, ls2 = overlaps.overlap_pairs(*args, other_type, maxgap, minoverlap)
                q2, s2 = qidx[lq2], sidx[ls2]
                qstr2, sstr2 = self._strand[q2], subj._strand[s2]
                ok2 = ((qstr2 == STAR) | (sstr2 == STAR) | (qstr2 == sstr2)) & ((qstr2 == MINUS) | (sstr2 == MINUS))
                q, s = np.concatenate((q[~minus], q2[ok2])), np.concatenate((s[~minus], s2[ok2]))
            out_q.append(q)
            out_s.append(s)
        if not out_q:
            return empty_int(), empty_int()
        q, s = np.concatenate(out_q), np.concatenate(out_s)
        o = np.lexsort((s, q))
        return q[o], s[o]

    def find_overlaps(
        self,
        subject: "GRanges | None" = None,
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",
        select: str = "all",
        ignore_strand: bool = False,
    ) -> Any:
        """As IRanges.find_overlaps, restricted to same seqlevel and compatible strand."""
        subj = self if subject is None else subject
        qi, si = self._overlap_pairs(subj, maxgap, minoverlap, type, ignore_strand)
        if select == "all":
            return Hits(qi, si, len(self), len(subj))
        return overlaps.select_hits(qi, si, len(self), select)

    def count_overlaps(
        self, subject: "GRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any", ignore_strand: bool = False
    ) -> IntArray:
        return self.find_overlaps(subject, maxgap, minoverlap, type, ignore_strand=ignore_strand).count_per_query()

    def overlaps_any(
        self, subject: "GRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any", ignore_strand: bool = False
    ) -> Any:
        return self.count_overlaps(subject, maxgap, minoverlap, type, ignore_strand) > 0

    def subset_by_overlaps(
        self,
        subject: "GRanges",
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",
        invert: bool = False,
        ignore_strand: bool = False,
    ) -> Self:
        mask = self.overlaps_any(subject, maxgap, minoverlap, type, ignore_strand)
        return self[~mask if invert else mask]

    # --- nearest -----------------------------------------------------------------
    # Same seqlevel and compatible strand required. For precede/follow a `-`
    # query looks toward lower coordinates; `*` queries are treated as `+`.

    def _nearest_family(
        self, kind: str, subject: "GRanges | None", select: str, allowed: tuple[str, str], ignore_strand: bool
    ) -> Any:
        if select not in allowed:
            raise ValueError(f"select must be one of {allowed}, got {select!r}")
        subj = self if subject is None else subject
        _, q_seq, s_seq = self._align(subj)
        qi, si = nearest.nearest_family(
            kind, self._ranges.start, self._ranges.end, q_seq, self._strand,
            subj._ranges.start, subj._ranges.end, s_seq, subj._strand,
            ignore_strand=ignore_strand, exclude_self=subject is None,
        )
        if select == "all":
            return Hits(qi, si, len(self), len(subj))
        return nearest.first_per_query(qi, si, len(self))

    def precede(self, subject: "GRanges | None" = None, select: str = "first", ignore_strand: bool = False) -> Any:
        return self._nearest_family("precede", subject, select, ("first", "all"), ignore_strand)

    def follow(self, subject: "GRanges | None" = None, select: str = "first", ignore_strand: bool = False) -> Any:
        return self._nearest_family("follow", subject, select, ("first", "all"), ignore_strand)

    def nearest(self, subject: "GRanges | None" = None, select: str = "arbitrary", ignore_strand: bool = False) -> Any:
        return self._nearest_family("nearest", subject, select, ("arbitrary", "all"), ignore_strand)

    def distance(self, other: "GRanges", ignore_strand: bool = False) -> IntArray:
        """Parallel distance; `NO_HIT` (-1, R `NA`) across seqlevels or incompatible strands."""
        if len(self) != len(other):
            raise ValueError("distance requires equal lengths")
        _, a_seq, b_seq = self._align(other)
        a, b = self._strand, other._strand
        d = nearest.distance(self._ranges.start, self._ranges.end, other._ranges.start, other._ranges.end)
        ok = a_seq == b_seq
        if not ignore_strand:
            ok &= (a == b) | (a == STAR) | (b == STAR)
        return np.where(ok, d, NO_HIT)

    # --- coverage ----------------------------------------------------------------

    def coverage(
        self, shift: int = 0, width: int | Mapping[str, int] | None = None, weight: Any = 1
    ) -> dict[str, Rle]:
        """Coverage per seqlevel (strand ignored) as `{seqname: Rle}`.

        Length defaults to the seqlength, else the largest shifted end.
        """
        w = np.asarray(weight)
        if w.ndim == 0:
            w = np.full(len(self), w)
        if len(w) != len(self):
            raise ValueError("weight must be scalar or one per range")
        out: dict[str, Rle] = {}
        for k, name in enumerate(self._seqinfo.seqnames):
            idx = np.flatnonzero(self._seq == k)
            if isinstance(width, Mapping):
                wd = width.get(name)
            else:
                wd = width
            if wd is None:
                wd = self._seqinfo.seqlengths[k]
            out[name] = coverage.coverage_kernel(
                self._ranges.start[idx], self._ranges.end[idx], shift, wd, w[idx]
            )
        return out
