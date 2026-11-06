"""GRanges: IRanges + seqnames + strand + seqinfo."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _group, _mcols, inter, intra, setops
from .iranges import IRanges, resolve_subscript
from .seqinfo import Seqinfo

STRANDS = ("+", "-", "*")
STRAND_CODE = {s: i for i, s in enumerate(STRANDS)}


def _recycled(values: Any, n: int, what: str) -> list[Any]:
    vals = [values] if isinstance(values, str) else list(values)
    if len(vals) == 1 and n != 1:
        return vals * n
    if len(vals) != n:
        raise ValueError(f"'{what}' has length {len(vals)}, expected 1 or {n}")
    return vals


class GRanges:
    """Genomic ranges. Per-range codes index into ``seqinfo``; strand is 0 +, 1 -, 2 *."""

    __slots__ = ("_seqcodes", "_strand", "_ranges", "_mcols", "_seqinfo")

    def __init__(
        self,
        seqnames: Sequence[str],
        ranges: IRanges | None = None,
        strand: Sequence[str] | None = None,
        seqinfo: Seqinfo | None = None,
        *,
        start: Any = None,
        end: Any = None,
        width: Any = None,
        names: Sequence[str | None] | None = None,
        mcols: Any = None,
    ) -> None:
        if ranges is None:
            ranges = IRanges(start=start, end=end, width=width, names=names)
        names_in = [seqnames] if isinstance(seqnames, str) else list(seqnames)
        n = max(len(names_in), len(ranges), len(strand) if strand is not None and not isinstance(strand, str) else 1)
        if len(names_in) == 0:
            n = len(ranges)
        sq = _recycled(names_in, n, "seqnames")
        if any(s is None for s in sq):
            raise ValueError("'seqnames' cannot contain NAs")
        if len(ranges) == 1 and n != 1:
            ranges = IRanges.concat(*([ranges] * n))
        if len(ranges) != n:
            raise ValueError(f"'ranges' must have length {n} or 1")
        st = ["*"] * n if strand is None else _recycled(strand, n, "strand")
        if any(s not in STRAND_CODE for s in st):
            raise ValueError("strand values must be in '+', '-', '*'")
        if seqinfo is None:
            seqinfo = Seqinfo(list(dict.fromkeys(sq)))
        elif not set(sq) <= set(seqinfo.seqnames):
            raise ValueError("'seqnames' contains sequence names with no entries in 'seqinfo'")
        lookup = {s: i for i, s in enumerate(seqinfo.seqnames)}
        self._init(
            np.array([lookup[s] for s in sq], dtype=np.int64),
            np.array([STRAND_CODE[s] for s in st], dtype=np.int8),
            ranges.with_mcols(None),
            _mcols.coerce(mcols, n),
            seqinfo,
        )

    def _init(self, codes: np.ndarray, strand: np.ndarray, ranges: IRanges, mcols: pl.DataFrame, seqinfo: Seqinfo) -> None:
        self._seqcodes = codes
        self._strand = strand
        self._ranges = ranges
        self._mcols = mcols
        self._seqinfo = seqinfo

    @classmethod
    def _make(cls, codes: np.ndarray, strand: np.ndarray, ranges: IRanges, mcols: pl.DataFrame | None, seqinfo: Seqinfo) -> Self:
        obj = object.__new__(cls)
        obj._init(codes.astype(np.int64), strand.astype(np.int8), ranges.with_mcols(None), mcols if mcols is not None else _mcols.empty(), seqinfo)
        return obj

    # accessors -----------------------------------------------------------
    @property
    def seqnames(self) -> np.ndarray:
        return np.array(self._seqinfo.seqnames, dtype=object)[self._seqcodes] if len(self) else np.zeros(0, dtype=object)

    @property
    def strand(self) -> np.ndarray:
        return np.array(STRANDS, dtype=object)[self._strand] if len(self) else np.zeros(0, dtype=object)

    @property
    def ranges(self) -> IRanges:
        return self._ranges

    @property
    def start(self) -> np.ndarray:
        return self._ranges.start

    @property
    def end(self) -> np.ndarray:
        return self._ranges.end

    @property
    def width(self) -> np.ndarray:
        return self._ranges.width

    @property
    def names(self) -> np.ndarray | None:
        return self._ranges.names

    @property
    def mcols(self) -> pl.DataFrame:
        return self._mcols

    @property
    def seqinfo(self) -> Seqinfo:
        return self._seqinfo

    @property
    def seqlevels(self) -> list[str]:
        return self._seqinfo.seqlevels

    @property
    def seqlengths(self) -> dict[str, int | None]:
        return self._seqinfo.seqlengths

    @property
    def is_circular(self) -> dict[str, bool | None]:
        return self._seqinfo.is_circular

    @property
    def genome(self) -> dict[str, str | None]:
        return self._seqinfo.genome

    def with_mcols(self, mcols: Any) -> Self:
        return self._make(self._seqcodes, self._strand, self._ranges, _mcols.coerce(mcols, len(self)), self._seqinfo)

    def __len__(self) -> int:
        return len(self._ranges)

    # subsetting / combining ----------------------------------------------
    def _take(self, idx: np.ndarray) -> Self:
        return self._make(
            self._seqcodes[idx], self._strand[idx], self._ranges._take(idx), _mcols.take(self._mcols, idx), self._seqinfo
        )

    def __getitem__(self, i: Any) -> Self:
        return self._take(resolve_subscript(i, len(self), self._ranges._names))

    @classmethod
    def concat(cls, *xs: Self) -> Self:
        if not xs:
            raise ValueError("concat needs at least one GRanges")
        si = xs[0]._seqinfo
        for x in xs[1:]:
            si = si.merge(x._seqinfo)
        lookup = {s: i for i, s in enumerate(si.seqnames)}
        codes = np.concatenate(
            [np.array([lookup[s] for s in x._seqinfo.seqnames], dtype=np.int64)[x._seqcodes] for x in xs]
        )
        return cls._make(
            codes,
            np.concatenate([x._strand for x in xs]),
            IRanges.concat(*(x._ranges for x in xs)),
            _mcols.concat([x._mcols for x in xs], [len(x) for x in xs]),
            si,
        )

    # intra-range ---------------------------------------------------------
    def _with_coords(self, start: np.ndarray, width: np.ndarray, use_names: bool = True) -> Self:
        names = self._ranges._names if use_names else None
        return self._make(self._seqcodes, self._strand, IRanges._make(start, width, names), self._mcols, self._seqinfo)

    def _reverse(self, ignore_strand: bool) -> np.ndarray:
        return np.zeros(len(self), bool) if ignore_strand else self._strand == 1

    def shift(self, shift: Any) -> Self:
        return self._with_coords(*intra.shift(self._ranges._start, self._ranges._width, intra.arg_vector(shift, len(self), "shift")))

    def narrow(self, start: Any = None, end: Any = None, width: Any = None) -> Self:
        r = self._ranges.narrow(start, end, width)
        return self._with_coords(r._start, r._width)

    def resize(self, width: Any, fix: Any = "start", use_names: bool = True, ignore_strand: bool = False) -> Self:
        n = len(self)
        return self._with_coords(
            *intra.resize(
                self._ranges._start, self._ranges._width, intra.arg_vector(width, n, "width"),
                intra.str_vector(fix, n, "fix"), self._reverse(ignore_strand),
            ),
            use_names,
        )

    def flank(self, width: Any, start: Any = True, both: bool = False, use_names: bool = True, ignore_strand: bool = False) -> Self:
        n = len(self)
        return self._with_coords(
            *intra.flank(
                self._ranges._start, self._ranges._width, intra.arg_vector(width, n, "width"),
                intra.flag_vector(start, n, "start"), both, self._reverse(ignore_strand),
            ),
            use_names,
        )

    def promoters(self, upstream: Any = 2000, downstream: Any = 200, use_names: bool = True) -> Self:
        n = len(self)
        return self._with_coords(
            *intra.promoters(
                self._ranges._start, self._ranges._width, intra.arg_vector(upstream, n, "upstream"),
                intra.arg_vector(downstream, n, "downstream"), self._reverse(False),
            ),
            use_names,
        )

    # inter-range ---------------------------------------------------------
    def _from_groups(self, out: _group.Collected, seqinfo: Seqinfo | None = None) -> Self:
        codes, strands, start, width = out.flat()
        mcols = None
        if out.has_members:
            mcols = pl.DataFrame({"revmap": pl.Series("revmap", [m.tolist() for m in out.members], dtype=pl.List(pl.Int64))})
        return self._make(codes, strands, IRanges._make(start, width), mcols, seqinfo or self._seqinfo)

    def _by_group(self, kernel: Any, ignore_strand: bool) -> Self:
        return self._from_groups(
            _group.by_group(self._seqcodes, self._strand, self._ranges._start, self._ranges._width, kernel, ignore_strand)
        )

    def range(self, ignore_strand: bool = False) -> Self:
        return self._by_group(lambda s, w: (*inter.range_of(s, w), None), ignore_strand)

    def reduce(
        self, drop_empty_ranges: bool = False, min_gapwidth: int = 1, with_revmap: bool = False, ignore_strand: bool = False
    ) -> Self:
        def kernel(s: np.ndarray, w: np.ndarray) -> _group.Result:
            rs, rw, members = inter.reduce_runs(s, w, int(min_gapwidth), drop_empty_ranges)
            return rs, rw, members if with_revmap else None

        return self._by_group(kernel, ignore_strand)

    def disjoin(self, with_revmap: bool = False, ignore_strand: bool = False) -> Self:
        def kernel(s: np.ndarray, w: np.ndarray) -> _group.Result:
            rs, rw, members = inter.disjoin(s, w)
            return rs, rw, members if with_revmap else None

        return self._by_group(kernel, ignore_strand)

    def gaps(self, start: int | None = None, end: int | None = None, ignore_strand: bool = False) -> Self:
        """Gaps per (seqlevel, strand) within ``[start, end]``; default ``[1, seqlength]``.

        Every seqlevel gets a group for each strand (just ``*`` when ignoring
        strand), so unused sequences yield whole-sequence gaps.
        """
        lo = 1 if start is None else start
        lengths = self._seqinfo._seqlengths
        strands = (2,) if ignore_strand else (0, 1, 2)
        all_keys = [c * 3 + st for c in range(len(self._seqinfo)) for st in strands]
        out = _group.Collected()
        k = _group.keys(self._seqcodes, self._strand, ignore_strand)
        for key in all_keys:
            idx = np.flatnonzero(k == key)
            hi = end if end is not None else lengths[key // 3]
            s, w = inter.gaps(self._ranges._start[idx], self._ranges._width[idx], lo, hi)
            out.add(key, s, w)
        return self._from_groups(out)

    # set operations ------------------------------------------------------
    def _set_op(self, y: Self, kernel: Any, ignore_strand: bool) -> Self:
        si = self._seqinfo.merge(y._seqinfo)
        lookup = {s: i for i, s in enumerate(si.seqnames)}

        def mapped(g: Self) -> np.ndarray:
            return np.array([lookup[s] for s in g._seqinfo.seqnames], dtype=np.int64)[g._seqcodes]

        xk = _group.keys(mapped(self), self._strand, ignore_strand)
        yk = _group.keys(mapped(y), y._strand, ignore_strand)
        out = _group.by_group_pair(xk, (self._ranges._start, self._ranges._width), yk, (y._ranges._start, y._ranges._width), kernel)
        return self._from_groups(out, si)

    def union(self, y: Self, ignore_strand: bool = False) -> Self:
        return self._set_op(y, setops.union, ignore_strand)

    def intersect(self, y: Self, ignore_strand: bool = False) -> Self:
        return self._set_op(y, setops.intersect, ignore_strand)

    def setdiff(self, y: Self, ignore_strand: bool = False) -> Self:
        return self._set_op(y, setops.setdiff, ignore_strand)

    # parallel set operations ---------------------------------------------
    def _pair_compat(self, y: Self, ignore_strand: bool) -> tuple[np.ndarray, Seqinfo, np.ndarray]:
        """(compatible mask, merged seqinfo, y's codes in merged seqinfo)."""
        if len(self) != len(y):
            raise ValueError("'x' and 'y' must have the same length")
        si = self._seqinfo.merge(y._seqinfo)
        lookup = {s: i for i, s in enumerate(si.seqnames)}
        xc = np.array([lookup[s] for s in self._seqinfo.seqnames], dtype=np.int64)[self._seqcodes] if len(self) else self._seqcodes
        yc = np.array([lookup[s] for s in y._seqinfo.seqnames], dtype=np.int64)[y._seqcodes] if len(self) else y._seqcodes
        ok = xc == yc
        if not ignore_strand:
            ok &= (self._strand == y._strand) | (self._strand == 2) | (y._strand == 2)
        return ok, si, xc

    def _pair_result(self, si: Seqinfo, xc: np.ndarray, strand: np.ndarray, start: np.ndarray, width: np.ndarray, mcols: pl.DataFrame | None) -> Self:
        return self._make(xc, strand, IRanges._make(start, width, self._ranges._names), mcols, si)

    def _resolved_strand(self, y: Self, ignore_strand: bool) -> np.ndarray:
        """Strand of a compatible pair: x's, or y's when x is `*` (unless ignoring strand)."""
        if ignore_strand:
            return self._strand
        return np.where(self._strand == 2, y._strand, self._strand)

    def punion(self, y: Self, fill_gap: bool = False, ignore_strand: bool = False) -> Self:
        ok, si, xc = self._pair_compat(y, ignore_strand)
        if not ok.all():
            raise ValueError("'x' and 'y' elements must have compatible 'seqnames' and 'strand' values")
        s, w = setops.punion(self._ranges._coords(), y._ranges._coords(), fill_gap)
        return self._pair_result(si, xc, self._strand, s, w, None)

    def pintersect(self, y: Self, ignore_strand: bool = False) -> Self:
        """Adds a logical ``hit`` mcol: compatible pair with a non-ambiguous intersection."""
        ok, si, xc = self._pair_compat(y, ignore_strand)
        xs, xe = self._ranges._start, self._ranges.end
        s, e = np.maximum(xs, y._ranges._start), np.minimum(xe, y._ranges.end)
        hit = ok & (e >= s - 1)
        start = np.where(hit, s, xs)
        width = np.where(hit, e - s + 1, 0)
        strand = np.where(hit, self._resolved_strand(y, ignore_strand), self._strand)
        hit_col = pl.Series("hit", hit.tolist(), dtype=pl.Boolean)
        mcols = self._mcols.with_columns(hit_col) if self._mcols.width else pl.DataFrame([hit_col])
        return self._pair_result(si, xc, strand, start, width, mcols)

    def psetdiff(self, y: Self, ignore_strand: bool = False) -> Self:
        ok, si, xc = self._pair_compat(y, ignore_strand)
        s, w = self._ranges._start.copy(), self._ranges._width.copy()
        if ok.any():
            xs = (s[ok], w[ok])
            ns, nw = setops.psetdiff(xs, (y._ranges._start[ok], y._ranges._width[ok]))
            s[ok], w[ok] = ns, nw
        strand = np.where(ok, self._resolved_strand(y, ignore_strand), self._strand)
        return self._pair_result(si, xc, strand, s, w, None)

    # seqinfo-aware -------------------------------------------------------
    def trim(self) -> Self:
        """Clamp ranges to ``[1, seqlength]`` on non-circular sequences of known length."""
        lengths = np.array([-1 if v is None else v for v in self._seqinfo._seqlengths], dtype=np.int64)
        circ = np.array([bool(v) for v in self._seqinfo._is_circular])
        if not len(self):
            return self
        L = lengths[self._seqcodes]
        act = (L >= 0) & ~circ[self._seqcodes]
        s, e = self._ranges._start, self._ranges.end
        ns = np.where(act, np.clip(s, 1, np.where(act, L + 1, 1)), s)
        ne = np.where(act, np.clip(e, 0, np.where(act, L, 0)), e)
        return self._make(
            self._seqcodes, self._strand, IRanges._make(ns, ne - ns + 1, self._ranges._names), self._mcols, self._seqinfo
        )

    def __repr__(self) -> str:
        rows = ", ".join(f"{a}:{b}-{c}({d})" for a, b, c, d in zip(self.seqnames, self.start, self.end, self.strand, strict=True))
        return f"GRanges({rows})"
