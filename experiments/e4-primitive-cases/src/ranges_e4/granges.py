"""GRanges: IRanges + seqnames + strand + seqinfo."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _group, _mcols, inter, intra, nearest, ordering, overlaps, setops
from . import coverage as _coverage
from .rle import Rle
from .hits import Hits
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

    # overlaps ------------------------------------------------------------
    def _aligned(self, subject: Self) -> tuple[np.ndarray, np.ndarray]:
        """Seqlevel codes of self and subject in their merged seqinfo (conflicts raise)."""
        si = self._seqinfo.merge(subject._seqinfo)
        lookup = {s: i for i, s in enumerate(si.seqnames)}

        def mapped(g: Self) -> np.ndarray:
            return np.array([lookup[s] for s in g._seqinfo.seqnames], dtype=np.int64)[g._seqcodes]

        return mapped(self), mapped(subject)

    def _overlap_pairs(
        self, subject: Self, maxgap: int, minoverlap: int, type: str, ignore_strand: bool
    ) -> tuple[np.ndarray, np.ndarray]:
        """Offset each seqlevel by 2^34 so the IRanges kernel never pairs different sequences."""
        qc, sc = self._aligned(subject)
        gap = 1 << 34
        qs, ss = self._ranges._start + qc * gap, subject._ranges._start + sc * gap
        qw, sw = self._ranges._width, subject._ranges._width
        if ignore_strand:
            combos = [(np.arange(len(self)), np.arange(len(subject)))]
        else:
            compat = [(0, 0), (0, 2), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
            combos = [(np.flatnonzero(self._strand == a), np.flatnonzero(subject._strand == b)) for a, b in compat]
        qs_out: list[np.ndarray] = []
        ss_out: list[np.ndarray] = []
        for qi, si in combos:
            q, s = overlaps.overlap_pairs(qs[qi], qw[qi], ss[si], sw[si], maxgap, minoverlap, type)
            qs_out.append(qi[q])
            ss_out.append(si[s])
        q, s = np.concatenate(qs_out), np.concatenate(ss_out)
        o = np.lexsort((s, q))
        return q[o], s[o]

    def find_overlaps(
        self,
        subject: Self | None = None,
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",
        select: str = "all",
        ignore_strand: bool = False,
    ) -> Hits | np.ndarray:
        subject = self if subject is None else subject
        q, s = self._overlap_pairs(subject, int(maxgap), int(minoverlap), type, ignore_strand)
        if select == "all":
            return Hits(q, s, len(self), len(subject))
        return overlaps.select_hits(q, s, len(self), select)

    def count_overlaps(
        self, subject: Self | None = None, maxgap: int = -1, minoverlap: int = 0, type: str = "any", ignore_strand: bool = False
    ) -> np.ndarray:
        subject = self if subject is None else subject
        q, _ = self._overlap_pairs(subject, int(maxgap), int(minoverlap), type, ignore_strand)
        return np.bincount(q, minlength=len(self)).astype(np.int64)

    def overlaps_any(
        self, subject: Self | None = None, maxgap: int = -1, minoverlap: int = 0, type: str = "any", ignore_strand: bool = False
    ) -> np.ndarray:
        return self.count_overlaps(subject, maxgap, minoverlap, type, ignore_strand) > 0

    def subset_by_overlaps(
        self,
        subject: Self | None = None,
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",
        invert: bool = False,
        ignore_strand: bool = False,
    ) -> Self:
        keep = self.overlaps_any(subject, maxgap, minoverlap, type, ignore_strand)
        return self[~keep if invert else keep]

    # nearest -------------------------------------------------------------
    def _pair_arrays(self, subject: Self, ignore_strand: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        qc, sc = self._aligned(subject)
        qstr = np.full(len(self), 2) if ignore_strand else self._strand
        sstr = np.full(len(subject), 2) if ignore_strand else subject._strand
        return qc, sc, qstr, sstr

    def _nearest_family(self, subject: Self, select: str, allowed: tuple[str, str], ignore_strand: bool, mode: str) -> Hits | np.ndarray:
        """mode: nearest | precede | follow. Eligible pairs share a seqlevel and have compatible strands.

        precede/follow direction is the query's strand, or the subject's when
        the query is `*` (`*` on both counts as `+`): on `-` the genome axis is flipped.
        """
        if select not in allowed:
            raise ValueError(f"'select' must be one of {allowed}")
        qc, sc, qstr, sstr = self._pair_arrays(subject, ignore_strand)
        qs, qe = self._ranges._start, self._ranges.end
        ss, se = subject._ranges._start, subject._ranges.end
        per_query: list[np.ndarray] = []
        for i in range(len(self)):
            eligible = (sc == qc[i]) & ((sstr == qstr[i]) | (sstr == 2) | (qstr[i] == 2))
            if mode == "nearest":
                dist = np.maximum(0, np.maximum(qs[i], ss) - np.minimum(qe[i], se) - 1)
            else:
                plus = (qstr[i] == 0) | ((qstr[i] == 2) & (sstr != 1))
                right = ss > qe[i]  # subject lies after the query on the genome axis
                left = se < qs[i]
                right_gap, left_gap = ss - qe[i] - 1, qs[i] - se - 1
                if mode == "precede":  # downstream: right on `+`, left on `-`
                    eligible &= np.where(plus, right, left)
                    dist = np.where(plus, right_gap, left_gap)
                else:  # upstream: left on `+`, right on `-`
                    eligible &= np.where(plus, left, right)
                    dist = np.where(plus, left_gap, right_gap)
            per_query.append(nearest.ties_at_minimum(dist, eligible))
        q, s = nearest.collect(per_query)
        if select == "all":
            return Hits(q, s, len(self), len(subject))
        return overlaps.select_hits(q, s, len(self), "last" if select == "last" else "first")

    def nearest(self, subject: Self, select: str = "arbitrary", ignore_strand: bool = False) -> Hits | np.ndarray:
        return self._nearest_family(subject, select, ("arbitrary", "all"), ignore_strand, "nearest")

    def precede(self, subject: Self, select: str = "first", ignore_strand: bool = False) -> Hits | np.ndarray:
        return self._nearest_family(subject, select, ("first", "all"), ignore_strand, "precede")

    def follow(self, subject: Self, select: str = "last", ignore_strand: bool = False) -> Hits | np.ndarray:
        return self._nearest_family(subject, select, ("last", "all"), ignore_strand, "follow")

    def distance(self, y: Self, ignore_strand: bool = False) -> np.ndarray:
        """Masked where the pair is on different seqlevels or has incompatible strands."""
        if len(self) != len(y) and 1 not in (len(self), len(y)):
            raise ValueError("'x' and 'y' must have the same length (or one of length 1)")
        n = max(len(self), len(y)) if len(self) and len(y) else 0
        qc, sc = self._aligned(y)
        bc = lambda a: np.broadcast_to(a, (n,))  # noqa: E731
        d = nearest.distance(bc(self._ranges._start), bc(self._ranges._width), bc(y._ranges._start), bc(y._ranges._width))
        ok = bc(qc) == bc(sc)
        if not ignore_strand:
            qs_, ys_ = bc(self._strand), bc(y._strand)
            ok = ok & ((qs_ == ys_) | (qs_ == 2) | (ys_ == 2))
        return np.ma.MaskedArray(d, mask=~ok)

    # coverage ------------------------------------------------------------
    def coverage(self, shift: Any = 0, width: Any = None, weight: Any = 1) -> dict[str, Rle]:
        """Per-seqlevel coverage (strand ignored), in seqlevels order.

        Length of each vector: ``width`` (scalar, per-seqlevel sequence, or dict
        by seqname) if given, else the seqlength, else the max end. ``weight``
        may be a scalar, a per-range vector, or the name of a numeric mcols column.
        """
        n, nlev = len(self), len(self._seqinfo)
        if isinstance(weight, str):
            if weight not in self._mcols.columns:
                raise ValueError(f"no mcols column named {weight!r}")
            weight = self._mcols[weight].to_numpy()
        w_arr = np.asarray(weight)
        shift_arr = intra.arg_vector(shift, n, "shift") if n else np.zeros(0, np.int64)
        if isinstance(width, dict):
            widths: list[int | None] = [width.get(s) for s in self._seqinfo.seqnames]
        elif width is None:
            widths = [None] * nlev
        else:
            vals = np.atleast_1d(np.asarray(width)).tolist()
            if len(vals) == 1:
                vals = vals * nlev
            if len(vals) != nlev:
                raise ValueError(f"'width' must have length 1 or {nlev}")
            widths = [None if v is None else int(v) for v in vals]
        out: dict[str, Rle] = {}
        for c, name in enumerate(self._seqinfo.seqnames):
            idx = np.flatnonzero(self._seqcodes == c)
            length = widths[c] if widths[c] is not None else self._seqinfo._seqlengths[c]
            wc = w_arr if w_arr.ndim == 0 or len(w_arr) == 1 else w_arr[idx]
            dense = _coverage.coverage_dense(
                self._ranges._start[idx], self._ranges._width[idx], shift_arr[idx], length, wc
            )
            out[name] = Rle.encode(dense)
        return out

    # ordering ------------------------------------------------------------
    def _keys(self) -> np.ndarray:
        """Sort key rows: seqlevel position, strand (+ - *), start, width."""
        if not len(self):
            return np.zeros((4, 0), np.int64)
        return np.vstack([self._seqcodes, self._strand.astype(np.int64), self._ranges._start, self._ranges._width])

    def order(self, decreasing: bool = False) -> np.ndarray:
        return ordering.order(self._keys(), decreasing)

    def sort(self, decreasing: bool = False) -> Self:
        return self._take(self.order(decreasing))

    def rank(self, *, ties_method: str) -> np.ndarray:
        return ordering.rank(self._keys(), ties_method)

    def duplicated(self, from_last: bool = False) -> np.ndarray:
        return ordering.duplicated(self._keys(), from_last)

    def unique(self) -> Self:
        return self[~self.duplicated()]

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
