"""GRanges: IRanges plus seqnames, strand and seqinfo."""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _args, _meta
from ._args import flags
from ._kernels import coverage as coverage_kernel
from ._kernels import inter, intra, overlaps, setops
from ._kernels import nearest as nearest_kernel
from ._kernels import order as order_kernel
from ._meta import Names
from .hits import Hits, overlap_result
from .iranges import IRanges
from .rle import Rle
from .seqinfo import Seqinfo

STRANDS = ("+", "-", "*")
PLUS, MINUS, STAR = 0, 1, 2


def _as_str_list(v: Any, what: str) -> list[str | None]:
    if isinstance(v, str):
        return [v]
    return [None if x is None else str(x) for x in list(v)]


def _parse_strand(v: Any) -> np.ndarray:
    vals = _as_str_list(v, "strand")
    out = np.empty(len(vals), dtype=np.int8)
    na = False
    for i, s in enumerate(vals):
        if s is None or s == "NA":
            out[i], na = STAR, True
        elif s in STRANDS:
            out[i] = STRANDS.index(s)
        else:
            raise ValueError("strand values must be in '+' '-' '*'")
    if na:
        warnings.warn("missing values in 'strand' converted to \"*\"", UserWarning, stacklevel=3)
    return out


class GRanges:
    """Genomic ranges. Coordinates are 1-based closed; seqinfo defines the sequence levels."""

    __slots__ = ("_seq", "_ranges", "_strand", "_names", "_mcols", "_seqinfo")

    def __init__(
        self,
        seqnames: Any,
        ranges: IRanges | None = None,
        *,
        start: Any = None,
        end: Any = None,
        width: Any = None,
        strand: Any = None,
        names: Iterable[str | None] | None = None,
        mcols: Any = None,
        seqinfo: Seqinfo | None = None,
        seqlengths: Any = None,
    ) -> None:
        if ranges is None:
            ranges = IRanges(start=start, end=end, width=width)
        elif not isinstance(ranges, IRanges):
            raise TypeError("'ranges' must be an IRanges")
        sn = _as_str_list(seqnames, "seqnames")
        if any(s is None for s in sn):
            raise ValueError("'seqnames' cannot contain NAs")
        st = _parse_strand("*" if strand is None else strand)
        n = max(len(sn), len(ranges), len(st))
        if 0 in (len(sn), len(ranges), len(st)):
            n = 0
        for what, k in (("seqnames", len(sn)), ("ranges", len(ranges)), ("strand", len(st))):
            if k != n and not (k == 1 and n > 0):
                raise ValueError(f"'{what}' must have the length of the object to construct ({n}) or length 1")
        if len(sn) == 1 and n != 1:
            sn = sn * n
        if seqinfo is None:
            levels = list(dict.fromkeys(s for s in sn if s is not None))
            seqinfo = Seqinfo(levels, seqlengths)
        elif seqlengths is not None:
            raise ValueError("supply either 'seqinfo' or 'seqlengths', not both")
        index = {s: i for i, s in enumerate(seqinfo.seqnames)}
        if any(s not in index for s in sn):
            raise ValueError("'seqnames' contains sequence names with no entries in 'seqinfo'")
        codes = np.array([index[s] for s in sn if s is not None], dtype=np.int32)
        r_names = ranges.names
        nm = _meta.check_names(names, n) if names is not None else r_names
        flat = IRanges._new(
            np.resize(ranges.start, n), np.resize(ranges.width, n), validate=False
        )
        self._build(codes, flat, np.resize(st, n), nm, _meta.check_mcols(mcols, n), seqinfo)

    def _build(
        self,
        seq: np.ndarray,
        ranges: IRanges,
        strand: np.ndarray,
        names: Names | None,
        mcols: pl.DataFrame | None,
        seqinfo: Seqinfo,
        *,
        check: bool = True,
    ) -> None:
        seq = np.ascontiguousarray(seq, dtype=np.int32)
        strand = np.ascontiguousarray(strand, dtype=np.int8)
        seq.flags.writeable = False
        strand.flags.writeable = False
        self._seq, self._ranges, self._strand = seq, ranges, strand
        self._names, self._mcols, self._seqinfo = names, mcols, seqinfo
        if check:
            self._warn_out_of_bound()

    @classmethod
    def _new(
        cls,
        seq: np.ndarray,
        start: np.ndarray,
        width: np.ndarray,
        strand: np.ndarray,
        names: Names | None,
        mcols: pl.DataFrame | None,
        seqinfo: Seqinfo,
        *,
        check: bool = True,
    ) -> Self:
        """Internal constructor from arrays (validates coordinates)."""
        obj = cls.__new__(cls)
        obj._build(seq, IRanges._new(start, width), strand, names, mcols, seqinfo, check=check)
        return obj

    # -- validity ------------------------------------------------------------
    def _oob_mask(self) -> np.ndarray:
        """Ranges outside a known, non-circular sequence length."""
        lens = np.array([-1 if (l is None or c) else l for l, c in zip(self._seqinfo.seqlengths, self._seqinfo.is_circular)] or [-1], dtype=np.int64)
        L = lens[self._seq]
        return (L >= 0) & ((self.start < 1) | (self.end > L))

    def _warn_out_of_bound(self) -> None:
        if len(self) and self._oob_mask().any():
            seqs = sorted({self._seqinfo.seqnames[i] for i in self._seq[self._oob_mask()]})
            warnings.warn(
                f"GRanges object contains {int(self._oob_mask().sum())} out-of-bound range(s) located on "
                f"sequence(s) {', '.join(seqs)}. Use trim() to trim these ranges.",
                UserWarning,
                stacklevel=4,
            )

    # -- accessors -----------------------------------------------------------
    @property
    def seqnames(self) -> np.ndarray:
        return np.array(self._seqinfo.seqnames, dtype=object)[self._seq] if len(self) else np.empty(0, dtype=object)

    @property
    def strand(self) -> np.ndarray:
        return np.array(STRANDS, dtype=object)[self._strand] if len(self) else np.empty(0, dtype=object)

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
    def ranges(self) -> IRanges:
        """Coordinates as IRanges, carrying names (not mcols)."""
        return self._ranges.set_names(self._names)

    @property
    def names(self) -> Names | None:
        return self._names

    @property
    def mcols(self) -> pl.DataFrame | None:
        return self._mcols

    @property
    def seqinfo(self) -> Seqinfo:
        return self._seqinfo

    @property
    def seqlevels(self) -> tuple[str, ...]:
        return self._seqinfo.seqnames

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
        rows = ", ".join(f"{a}:{s}-{e}({t})" for a, s, e, t in zip(self.seqnames[:4], self.start[:4], self.end[:4], self.strand[:4]))
        return f"GRanges(n={n}: {rows}{', ...' if n > 4 else ''})"

    # -- setters (copies) ----------------------------------------------------
    def _copy(self, **kw: Any) -> Self:
        obj = self.__class__.__new__(self.__class__)
        obj._build(
            kw.get("seq", self._seq), kw.get("ranges", self._ranges), kw.get("strand", self._strand),
            kw.get("names", self._names), kw.get("mcols", self._mcols), kw.get("seqinfo", self._seqinfo),
            check=kw.get("check", False),
        )
        return obj

    def set_names(self, names: Iterable[str | None] | None) -> Self:
        return self._copy(names=_meta.check_names(names, len(self)))

    def set_mcols(self, mcols: Any) -> Self:
        return self._copy(mcols=_meta.check_mcols(mcols, len(self)))

    def set_seqinfo(self, seqinfo: Seqinfo) -> Self:
        """Replace seqinfo; the new one must contain every sequence in use (same names)."""
        if not isinstance(seqinfo, Seqinfo):
            raise TypeError("'seqinfo' must be a Seqinfo")
        index = {s: i for i, s in enumerate(seqinfo.seqnames)}
        used = {self._seqinfo.seqnames[i] for i in self._seq}
        if not used <= set(index):
            raise ValueError("the supplied seqinfo lacks sequences in use")
        remap = np.array([index.get(s, -1) for s in self._seqinfo.seqnames], dtype=np.int32)
        return self._copy(seq=remap[self._seq], seqinfo=seqinfo, check=True)

    def _with_seqinfo_fields(self, what: str, value: Any) -> Self:
        si = self._seqinfo
        if isinstance(value, dict):
            if list(value) != list(si.seqnames):
                raise ValueError(f"the names on the supplied '{what}' vector must be identical to the seqnames")
            value = [value[s] for s in si.seqnames]
        fields = {"seqlengths": si.seqlengths, "is_circular": si.is_circular, "genome": si.genome}
        fields[what] = value
        return self.set_seqinfo(Seqinfo(si.seqnames, **fields))

    def set_seqlengths(self, value: Any) -> Self:
        return self._with_seqinfo_fields("seqlengths", value)

    def set_is_circular(self, value: Any) -> Self:
        return self._with_seqinfo_fields("is_circular", value)

    def set_genome(self, value: Any) -> Self:
        return self._with_seqinfo_fields("genome", value)

    # -- subsetting / combining ----------------------------------------------
    def _take(self, pos: np.ndarray) -> Self:
        return self._copy(
            seq=self._seq[pos],
            ranges=self._ranges._take(pos),
            strand=self._strand[pos],
            names=_meta.take_names(self._names, pos),
            mcols=_meta.take_mcols(self._mcols, pos),
        )

    def __getitem__(self, i: Any) -> Self:
        return self._take(_meta.resolve_index(i, len(self), self._names))

    @classmethod
    def concat(cls, parts: Sequence[GRanges]) -> Self:
        """Combine; seqinfo is merged (union of levels, NA filled, conflicts error)."""
        parts = list(parts)
        if not parts or not all(isinstance(p, GRanges) for p in parts):
            raise TypeError("concat needs a non-empty sequence of GRanges")
        si = Seqinfo.merge_all([p._seqinfo for p in parts])
        index = {s: i for i, s in enumerate(si.seqnames)}
        seqs = [np.array([index[s] for s in p._seqinfo.seqnames], dtype=np.int32)[p._seq] for p in parts]
        return cls._new(
            np.concatenate(seqs),
            np.concatenate([p.start for p in parts]),
            np.concatenate([p.width for p in parts]),
            np.concatenate([p._strand for p in parts]),
            _meta.concat_names([(p._names, len(p)) for p in parts]),
            _meta.concat_mcols([(p._mcols, len(p)) for p in parts]),
            si,
        )

    # -- trim ----------------------------------------------------------------
    @flags("use_names")
    def trim(self, use_names: bool = True) -> Self:
        """Clip ranges to [1, seqlength] on sequences with a known, non-circular length."""
        lens = np.array([-1 if (l is None or c) else l for l, c in zip(self._seqinfo.seqlengths, self._seqinfo.is_circular)] or [-1], dtype=np.int64)
        L = lens[self._seq]
        known = L >= 0
        new_start = np.where(known, np.minimum(np.maximum(self.start, 1), L + 1), self.start)
        new_end = np.where(known, np.maximum(np.minimum(self.end, L), new_start - 1), self.end)
        return self._new(
            self._seq, new_start, new_end - new_start + 1, self._strand,
            self._names if use_names else None, self._mcols, self._seqinfo, check=False,
        )

    # -- intra-range ---------------------------------------------------------
    def _neg(self, ignore_strand: bool) -> np.ndarray:
        """Ranges to be treated as minus-strand (their 5' end is the right end)."""
        return np.zeros(len(self), dtype=bool) if ignore_strand else self._strand == MINUS

    def _moved(self, ranges: IRanges, use_names: bool) -> Self:
        """Same sequences/strands/mcols with new coordinates; warns if out of bounds."""
        return self._copy(ranges=ranges.set_names(None).set_mcols(None), names=self._names if use_names else None, check=True)

    @flags("use_names")
    def shift(self, shift: Any = 0, use_names: bool = True) -> Self:
        """Shift along the genome; strand is ignored (R)."""
        return self._moved(self._ranges.shift(shift), use_names)

    @flags("use_names")
    def narrow(self, start: Any = None, end: Any = None, width: Any = None, use_names: bool = True) -> Self:
        """Narrow relative to each range's left end; strand is ignored (R)."""
        return self._moved(self._ranges.narrow(start, end, width), use_names)

    @flags("use_names", "ignore_strand")
    def resize(self, width: Any, fix: Any = "start", use_names: bool = True, ignore_strand: bool = False) -> Self:
        """Resize keeping the `fix` end; "start"/"end" refer to the transcript sense, so they swap on "-"."""
        r = self._ranges
        new_w = r._int_arg(width, "width", longer="truncate")
        if (new_w < 0).any():
            raise ValueError("'width' values must be non-negative")
        codes = r._fix_codes(fix)
        neg = self._neg(ignore_strand)
        swapped = np.where(codes == intra.FIX_START, intra.FIX_END, np.where(codes == intra.FIX_END, intra.FIX_START, codes))
        s, w = intra.resize(r.start, r.width, new_w, np.where(neg, swapped, codes))
        return self._moved(IRanges._new(s, w), use_names)

    @flags("both", "use_names", "ignore_strand")
    def flank(
        self, width: Any, start: Any = True, both: bool = False, use_names: bool = True, ignore_strand: bool = False
    ) -> Self:
        """Flanks upstream (start=True) or downstream of each range, in the strand's sense."""
        r = self._ranges
        by = r._int_arg(width, "width")
        flags = np.asarray(start)
        if flags.dtype.kind != "b":
            raise ValueError("'start' must be logical without NA's")
        at_start = _meta.recycle_arg(np.atleast_1d(flags), len(self), "start", longer="truncate")
        s, w = intra.flank(r.start, r.width, by, at_start ^ self._neg(ignore_strand), bool(both))
        return self._moved(IRanges._new(s, w), use_names)

    @flags("use_names")
    def promoters(self, upstream: Any = 2000, downstream: Any = 200, use_names: bool = True) -> Self:
        """Promoter windows around the start of each range, strand-aware."""
        r = self._ranges
        up = r._int_arg(upstream, "upstream", longer="truncate")
        down = r._int_arg(downstream, "downstream", longer="truncate")
        s, w = intra.promoters(r.start, r.width, up, down, self._strand == MINUS)
        return self._moved(IRanges._new(s, w), use_names)

    # -- inter-range -----------------------------------------------------------
    def _gid(self, ignore_strand: bool) -> np.ndarray:
        """Group id = seqlevel code * 3 + strand code (strand forced to "*" if ignored)."""
        st = np.full(len(self), STAR, dtype=np.int64) if ignore_strand else self._strand.astype(np.int64)
        return self._seq.astype(np.int64) * 3 + st

    def _from_groups(self, gid: np.ndarray, start: np.ndarray, width: np.ndarray, mcols: pl.DataFrame | None = None,
                     seqinfo: Seqinfo | None = None) -> Self:
        return self._new(
            (gid // 3).astype(np.int32), start, width, (gid % 3).astype(np.int8), None, mcols,
            seqinfo if seqinfo is not None else self._seqinfo, check=False,
        )

    @flags("ignore_strand")
    def range(self, *others: GRanges, ignore_strand: bool = False) -> Self:
        """Smallest range per (seqname, strand) group present."""
        x = self.concat([self, *others]) if others else self
        red = inter.range_(x._gid(ignore_strand), x.start, x.width)
        return x._from_groups(red.group, red.start, red.width)

    @flags("drop_empty_ranges", "with_revmap", "ignore_strand")
    def reduce(
        self, drop_empty_ranges: bool = False, min_gapwidth: int = 1, with_revmap: bool = False,
        ignore_strand: bool = False,
    ) -> Self:
        gap = _args.scalar_int(min_gapwidth, "min.gapwidth")
        assert gap is not None
        if gap < 0:
            raise ValueError("'min.gapwidth' must be non-negative")
        red = inter.reduce(self._gid(ignore_strand), self.start, self.width, gap, drop_empty_ranges)
        mc = _meta.revmap_frame(red.members) if with_revmap else None
        return self._from_groups(red.group, red.start, red.width, mc)

    @flags("with_revmap", "ignore_strand")
    def disjoin(self, with_revmap: bool = False, ignore_strand: bool = False) -> Self:
        red = inter.disjoin(self._gid(ignore_strand), self.start, self.width, with_revmap)
        mc = _meta.revmap_frame(red.members) if with_revmap else None
        return self._from_groups(red.group, red.start, red.width, mc)

    def _per_seqlevel(self, v: Any, what: str) -> list[int | None]:
        """Scalar / dict by seqlevel / sequence of per-seqlevel values -> list (None = NA)."""
        levels = self._seqinfo.seqnames
        if isinstance(v, dict):
            if list(v) != list(levels):
                raise ValueError(f"the names of '{what}' must be identical to the seqlevels")
            vals = list(v.values())
        elif v is None or np.isscalar(v):
            vals = [v] * len(levels)
        else:
            vals = np.ma.getdata(v).tolist() if not isinstance(v, list) else v
            if np.ma.isMaskedArray(v):
                vals = [None if m else x for x, m in zip(np.ma.getdata(v).tolist(), np.ma.getmaskarray(v))]
            if not vals:
                raise ValueError(f"'{what}' has no elements")
            if len(levels) % len(vals):
                warnings.warn(f"'{what}' length is not a multiple of the number of seqlevels", UserWarning, stacklevel=3)
            vals = [vals[i % len(vals)] for i in range(len(levels))]  # R recycles
        if any(isinstance(x, (bool, np.bool_, str)) for x in vals):
            raise ValueError(f"'{what}' must be a vector of integers")
        return [None if (x is None or (isinstance(x, float) and x != x)) else int(x) for x in vals]

    @flags("ignore_strand")
    def gaps(self, start: Any = 1, end: Any = "seqlengths", ignore_strand: bool = False) -> Self:
        """Uncovered stretches per (seqlevel, strand) within [start, end].

        `end="seqlengths"` (default) uses the sequence lengths. A bound that is None (NA) falls
        back, per group, to the extent of its non-empty ranges; a group with neither ranges
        nor bounds yields nothing. Groups without ranges still report whole-sequence gaps.
        """
        if isinstance(end, str):
            if end != "seqlengths":
                raise ValueError("'end' must be a number, None or \"seqlengths\"")
            end = list(self._seqinfo.seqlengths)
        lo_l, hi_l = self._per_seqlevel(start, "start"), self._per_seqlevel(end, "end")
        strands = [STAR] if ignore_strand else [PLUS, MINUS, STAR]
        n_lev = len(self._seqinfo)
        ngroups = n_lev * 3
        gid = self._gid(ignore_strand)
        ne = self.width > 0
        lo = np.ones(ngroups, dtype=np.int64)
        hi = np.zeros(ngroups, dtype=np.int64)  # lo > hi: no gaps
        active = np.zeros(ngroups, dtype=bool)
        for lev in range(n_lev):
            for st in strands:
                g = lev * 3 + st
                sel = ne & (gid == g)
                a, b = lo_l[lev], hi_l[lev]
                if a is None:
                    a = int(self.start[sel].min()) if sel.any() else None
                if b is None:
                    b = int(self.end[sel].max()) if sel.any() else None
                if a is not None and b is not None:
                    lo[g], hi[g], active[g] = a, b, True
        g_out, s_out, w_out = inter.gaps(gid, self.start, self.width, lo, hi)
        keep = active[g_out]
        return self._from_groups(g_out[keep], s_out[keep], w_out[keep])

    # -- set operations ------------------------------------------------------
    def _merged_with(self, other: GRanges) -> tuple[Seqinfo, np.ndarray, np.ndarray]:
        """Merged seqinfo and both objects' seqlevel codes remapped into it."""
        si = self._seqinfo.merge(other._seqinfo)
        index = {s: i for i, s in enumerate(si.seqnames)}
        a = np.array([index[s] for s in self._seqinfo.seqnames], dtype=np.int64)
        b = np.array([index[s] for s in other._seqinfo.seqnames], dtype=np.int64)
        return si, a[self._seq], b[other._seq]

    def _set_op(self, other: GRanges, fn: Any, ignore_strand: bool) -> Self:
        si, sa, sb = self._merged_with(other)
        st_a = np.full(len(self), STAR, dtype=np.int64) if ignore_strand else self._strand.astype(np.int64)
        st_b = np.full(len(other), STAR, dtype=np.int64) if ignore_strand else other._strand.astype(np.int64)
        g, s, w = fn(sa * 3 + st_a, self.start, self.width, sb * 3 + st_b, other.start, other.width)
        return self._from_groups(g, s, w, None, si)

    @flags("ignore_strand")
    def union(self, other: GRanges, ignore_strand: bool = False) -> Self:
        return self._set_op(other, setops.union, ignore_strand)

    @flags("ignore_strand")
    def intersect(self, other: GRanges, ignore_strand: bool = False) -> Self:
        return self._set_op(other, setops.intersect, ignore_strand)

    @flags("ignore_strand")
    def setdiff(self, other: GRanges, ignore_strand: bool = False) -> Self:
        return self._set_op(other, setops.setdiff, ignore_strand)

    # -- parallel set operations -----------------------------------------------
    def _pair_with(
        self, other: GRanges, *, recycle_other: bool, ignore_strand: bool, strict_strand: bool = False
    ) -> tuple[Seqinfo, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Merged seqinfo, other's (recycled) starts/widths/strands, and the pair compatibility mask."""
        n = len(self)
        if len(other) != n and not (recycle_other and len(other) == 1):
            raise ValueError("'y' must have the length of 'x' or length 1" if recycle_other else "'x' and 'y' must have the same length")
        si, sa, sb = self._merged_with(other)
        rep = np.zeros(n, dtype=np.int64) if len(other) == 1 and n != 1 else np.arange(n)
        sb, st_b = sb[rep], other._strand[rep]
        same_seq = sa == sb
        if ignore_strand:
            ok = same_seq
        elif strict_strand:
            ok = same_seq & (self._strand == st_b)
        else:
            ok = same_seq & ((self._strand == st_b) | (self._strand == STAR) | (st_b == STAR))
        return si, other.start[rep], other.width[rep], st_b, ok

    def _merged_strand(self, st_b: np.ndarray, ignore_strand: bool) -> np.ndarray:
        """Pair strand: x's, unless x is "*" (then y's) -- whether or not the pair is compatible."""
        if ignore_strand:
            return self._strand
        return np.where(self._strand == STAR, st_b, self._strand).astype(np.int8)

    @flags("fill_gap", "ignore_strand")
    def punion(self, other: GRanges, fill_gap: bool = False, ignore_strand: bool = False) -> Self:
        si, ys, yw, st_b, ok = self._pair_with(other, recycle_other=False, ignore_strand=ignore_strand)
        if not ok.all():
            raise ValueError("'x' and 'y' elements must have compatible 'seqnames' and 'strand' values")
        s, w = setops.punion(self.start, self.width, ys, yw, bool(fill_gap))
        return self._new(self._seq_in(si), s, w, self._strand, self._names, None, si, check=True)

    def _restricted(self, merged: Seqinfo) -> Seqinfo:
        """`merged` limited to this object's own seqlevels, in its order (what pintersect returns)."""
        idx = [merged.index(s) for s in self._seqinfo.seqnames]
        return Seqinfo(
            self._seqinfo.seqnames, [merged.seqlengths[i] for i in idx],
            [merged.is_circular[i] for i in idx], [merged.genome[i] for i in idx],
        )

    def _seq_in(self, si: Seqinfo) -> np.ndarray:
        """This object's seqlevel codes expressed in the (merged) seqinfo `si`."""
        index = {s: i for i, s in enumerate(si.seqnames)}
        return np.array([index[s] for s in self._seqinfo.seqnames], dtype=np.int32)[self._seq]

    @flags("drop_nohit_ranges", "ignore_strand", "strict_strand")
    def pintersect(
        self, other: GRanges, drop_nohit_ranges: bool = False, ignore_strand: bool = False, strict_strand: bool = False
    ) -> Self:
        """Pairwise intersection with a `hit` mcols column (False: incompatible or disjoint pair)."""
        si, ys, yw, st_b, ok = self._pair_with(
            other, recycle_other=True, ignore_strand=ignore_strand, strict_strand=strict_strand
        )
        s = np.maximum(self.start, ys)
        w = np.minimum(self.end, ys + yw - 1) - s + 1
        hit = ok & (w >= 0)
        start = np.where(hit, s, self.start)
        width = np.where(hit, w, 0)
        strand = np.where(hit, self._merged_strand(st_b, ignore_strand), self._strand).astype(np.int8)
        names, mcols = self._names, self._mcols
        out_si = self._restricted(si)
        if drop_nohit_ranges:
            keep = np.flatnonzero(hit)
            return self._new(
                self._seq[keep], start[keep], width[keep], strand[keep],
                _meta.take_names(names, keep), _meta.take_mcols(mcols, keep), out_si, check=True,
            )
        hit_col = pl.DataFrame([pl.Series("hit", hit.tolist(), dtype=pl.Boolean)])
        mc = hit_col if mcols is None else pl.concat([mcols, hit_col], how="horizontal")
        return self._new(self._seq, start, width, strand, names, mc, out_si, check=True)

    @flags("ignore_strand")
    def psetdiff(self, other: GRanges, ignore_strand: bool = False) -> Self:
        si, ys, yw, st_b, ok = self._pair_with(other, recycle_other=False, ignore_strand=ignore_strand)
        start, width = self.start.copy(), self.width.copy()
        idx = np.flatnonzero(ok)
        if len(idx):
            s, w = setops.psetdiff(start[idx], width[idx], ys[idx], yw[idx])
            start[idx], width[idx] = s, w
        strand = self._merged_strand(st_b, ignore_strand)
        return self._new(self._seq_in(si), start, width, strand, self._names, None, si, check=True)

    # -- overlaps ----------------------------------------------------------------
    _LIFT = 1 << 33  # separates seqlevels on one axis, as in the set-op kernels

    def _union_codes(self, other: GRanges) -> tuple[np.ndarray, np.ndarray]:
        """Seqlevel codes of both objects over the union of their seqlevels (no merge checks)."""
        levels = list(dict.fromkeys([*self._seqinfo.seqnames, *other._seqinfo.seqnames]))
        index = {s: i for i, s in enumerate(levels)}
        a = np.array([index[s] for s in self._seqinfo.seqnames], dtype=np.int64)[self._seq]
        b = np.array([index[s] for s in other._seqinfo.seqnames], dtype=np.int64)[other._seq]
        return a, b

    def _strand_ok(self, other_strand: np.ndarray, q: np.ndarray, s: np.ndarray, ignore_strand: bool) -> np.ndarray:
        if ignore_strand:
            return np.ones(len(q), dtype=bool)
        a, b = self._strand[q], other_strand[s]
        return (a == b) | (a == STAR) | (b == STAR)

    def _pairs(
        self, subject: GRanges, maxgap: int, minoverlap: int, type_: str, ignore_strand: bool
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Overlap pairs plus, per query, whether R's engine walked an NCList built on the query.

        R runs the engine per seqlevel and builds the NCList on the shorter side of each group
        (ties -> subject), which fixes hit order and `select="arbitrary"`.
        """
        si, qa, sb = self._merged_with(subject)
        n_lev = len(si)
        nq = np.bincount(qa, minlength=n_lev)
        ns = np.bincount(sb, minlength=n_lev)
        q, s = overlaps.find(
            self.start + qa * self._LIFT, self.width, subject.start + sb * self._LIFT, subject.width,
            maxgap, minoverlap, type_, ncl_order=(nq >= ns)[sb],
        )
        ok = self._strand_ok(subject._strand, q, s, ignore_strand)
        return q[ok], s[ok], (nq < ns)[qa]

    @flags("ignore_strand", "drop_self", "drop_redundant")
    def find_overlaps(
        self,
        subject: GRanges | None = None,
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",  # noqa: A002 - R's argument name
        select: str = "all",
        ignore_strand: bool = False,
        drop_self: bool = False,
        drop_redundant: bool = False,
    ) -> Hits | np.ma.MaskedArray:
        """Overlaps with `subject` (None = self); same seqname and compatible strand required."""
        type, select, gap, minov = overlaps.normalize_args(type, select, maxgap, minoverlap)
        subj = self if subject is None else subject
        q, s, shorter = self._pairs(subj, gap, minov, type, ignore_strand)
        return overlap_result(
            q, s, len(self), len(subj), select=select, self_hits=subject is None,
            drop_self=drop_self, drop_redundant=drop_redundant, query_is_shorter=shorter,
        )

    @flags("ignore_strand")
    def count_overlaps(
        self, subject: GRanges | None = None, maxgap: int = -1, minoverlap: int = 0,
        type: str = "any", ignore_strand: bool = False,  # noqa: A002
    ) -> np.ndarray:
        hits = self.find_overlaps(subject, maxgap, minoverlap, type, ignore_strand=ignore_strand)
        assert isinstance(hits, Hits)
        return hits.count_by_query()

    @flags("ignore_strand")
    def overlaps_any(
        self, subject: GRanges | None = None, maxgap: int = -1, minoverlap: int = 0,
        type: str = "any", ignore_strand: bool = False,  # noqa: A002
    ) -> np.ndarray:
        return self.count_overlaps(subject, maxgap, minoverlap, type, ignore_strand) > 0

    @flags("invert", "ignore_strand")
    def subset_by_overlaps(
        self, subject: GRanges, maxgap: int = -1, minoverlap: int = 0, type: str = "any",  # noqa: A002
        invert: bool = False, ignore_strand: bool = False,
    ) -> Self:
        hit = self.overlaps_any(subject, maxgap, minoverlap, type, ignore_strand)
        return self[~hit if invert else hit]

    # -- nearest / precede / follow / distance -------------------------------------
    def _geometry_hits(
        self, qpos: np.ndarray, spool: np.ndarray, qa: np.ndarray, sb: np.ndarray, subject: GRanges,
        precede_geometry: bool,
    ) -> tuple[np.ndarray, np.ndarray]:
        """(query, subject) pairs, all ties (ascending subject index): first subject starting after each query ends (precede geometry)
        or last subject ending before it starts (follow geometry), on the same seqlevel."""
        empty = np.empty(0, dtype=np.int64)
        if len(qpos) == 0 or len(spool) == 0:
            return empty, empty
        L = self._LIFT
        if precede_geometry:
            qi, si = nearest_kernel.precede(self.end[qpos] + qa[qpos] * L, subject.start[spool] + sb[spool] * L, True)
        else:
            qi, si = nearest_kernel.follow(self.start[qpos] + qa[qpos] * L, subject.end[spool] + sb[spool] * L, True)
        q, s = qpos[qi], spool[si]
        same = qa[q] == sb[s]
        # R quirk: it offsets each seqlevel by k*maxend and uses k*maxend as sentinels, so a subject
        # coordinate of 0 (an empty range at 1 for "follow", start 0 for "precede") on any seqlevel
        # but the first collides with a sentinel and the whole hit is dropped.
        coord = subject.start[s] if precede_geometry else subject.end[s]
        same &= ~((coord == 0) & (sb[s] >= 1))
        return q[same], s[same]

    def _precede_follow(
        self, subject: GRanges, where: str, ignore_strand: bool
    ) -> tuple[np.ndarray, np.ndarray]:
        """R's .GenomicRanges_findPrecedeFollow: all tied hits, sorted by (query, subject)."""
        empty = np.empty(0, dtype=np.int64)
        if len(self) == 0 or len(subject) == 0:
            return empty, empty
        qa, sb = self._union_codes(subject)
        qst = np.full(len(self), PLUS, dtype=np.int8) if ignore_strand else self._strand
        sst = np.full(len(subject), PLUS, dtype=np.int8) if ignore_strand else subject._strand
        plus_geom_precede = where == "precede"
        s_plus = np.flatnonzero(sst != MINUS)
        s_minus_incl = np.flatnonzero(sst != PLUS)
        s_minus_only = np.flatnonzero(sst == MINUS)
        q_plus, q_minus, q_star = (np.flatnonzero(qst == c) for c in (PLUS, MINUS, STAR))
        pq, ps = self._geometry_hits(q_plus, s_plus, qa, sb, subject, plus_geom_precede)
        mq, ms = self._geometry_hits(q_minus, s_minus_incl, qa, sb, subject, not plus_geom_precede)
        spq, sps = self._geometry_hits(q_star, s_plus, qa, sb, subject, plus_geom_precede)
        smq, sms = self._geometry_hits(q_star, s_minus_only, qa, sb, subject, not plus_geom_precede)
        both = np.intersect1d(spq, smq)
        sq, ss = np.concatenate([spq, smq]), np.concatenate([sps, sms])
        if len(both):
            gap = np.maximum(np.maximum(self.start[sq], subject.start[ss]) - np.minimum(self.end[sq], subject.end[ss]) - 1, 0)
            best = {int(b): int(gap[sq == b].min()) for b in both}
            keep = np.array([int(a) not in best or int(g) == best[int(a)] for a, g in zip(sq, gap)], dtype=bool)
            sq, ss = sq[keep], ss[keep]
        q_all, s_all = np.concatenate([pq, mq, sq]), np.concatenate([ps, ms, ss])
        order = np.lexsort((s_all, q_all))  # R's Hits come back sorted by (query, subject)
        return q_all[order], s_all[order]

    def _pick_one(self, q: np.ndarray, s: np.ndarray, last: bool) -> np.ma.MaskedArray:
        out = np.zeros(len(self), dtype=np.int64)
        has = np.zeros(len(self), dtype=bool)
        if len(q):
            if last:
                out[q] = s  # last write wins
            else:
                out[q[::-1]] = s[::-1]
            has[q] = True
        return np.ma.masked_array(out, mask=~has)

    @flags("ignore_strand")
    def precede(self, subject: GRanges | None = None, select: str = "first", ignore_strand: bool = False) -> Hits | np.ma.MaskedArray:
        select = _args.match_arg(select, ("first", "all"))
        subj = self if subject is None else subject
        q, s = self._precede_follow(subj, "precede", ignore_strand)
        return Hits(q, s, len(self), len(subj)) if select == "all" else self._pick_one(q, s, last=False)

    @flags("ignore_strand")
    def follow(self, subject: GRanges | None = None, select: str = "last", ignore_strand: bool = False) -> Hits | np.ma.MaskedArray:
        select = _args.match_arg(select, ("last", "all"))
        subj = self if subject is None else subject
        q, s = self._precede_follow(subj, "follow", ignore_strand)
        return Hits(q, s, len(self), len(subj)) if select == "all" else self._pick_one(q, s, last=True)

    @flags("ignore_strand")
    def nearest(
        self, subject: GRanges | None = None, select: str = "arbitrary", ignore_strand: bool = False
    ) -> Hits | np.ma.MaskedArray:
        """Overlap (maxgap 0) first, else the nearer of precede/follow; ties favour the lower index."""
        select = _args.match_arg(select, ("arbitrary", "all"))
        subj = self if subject is None else subject
        all_ = select == "all"
        ol = self.find_overlaps(subject, maxgap=0, select=select, ignore_strand=ignore_strand, drop_self=subject is None)
        olv = ol.select("first") if isinstance(ol, Hits) else ol
        miss = np.flatnonzero(np.ma.getmaskarray(olv))
        if len(miss) == 0:
            return ol
        sub = self._take(miss)
        pq, ps = sub._precede_follow(subj, "precede", ignore_strand)
        fq, fs = sub._precede_follow(subj, "follow", ignore_strand)
        if len(pq) == 0 and len(fq) == 0 and bool(np.ma.getmaskarray(olv).all()):
            return Hits(np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64), len(self), len(subj)) if all_ else olv
        p = sub._pick_one(pq, ps, last=False)  # first precede / last follow per query
        f = sub._pick_one(fq, fs, last=True)
        pm, fm = np.ma.getmaskarray(p), np.ma.getmaskarray(f)

        def dist(idx: np.ma.MaskedArray, m: np.ndarray) -> np.ndarray:
            i = np.where(m, 0, np.ma.getdata(idx))
            return np.where(m, -1, np.maximum(np.maximum(sub.start, subj.start[i]) - np.minimum(sub.end, subj.end[i]) - 1, 0))

        pd, fd = dist(p, pm), dist(f, fm)
        both = ~pm & ~fm
        pdata, fdata = np.ma.getdata(p), np.ma.getdata(f)
        pnearest = np.where(both, np.where(pd == fd, pdata < fdata, pd < fd), fm)
        if not all_:
            res = np.ma.getdata(olv).copy()
            mask = np.ma.getmaskarray(olv).copy()
            res[miss] = np.where(pnearest, pdata, fdata)
            mask[miss] = np.where(pnearest, pm, fm)
            return np.ma.masked_array(res, mask=mask)
        assert isinstance(ol, Hits)
        tie = both & (pd == fd)
        take_p = pnearest | tie
        take_f = ~pnearest & ~tie
        kp = take_p[pq]
        kf = take_f[fq]
        qs = np.concatenate([ol.query_hits, miss[pq[kp]], miss[fq[kf]]])
        ss = np.concatenate([ol.subject_hits, ps[kp], fs[kf]])
        order = np.lexsort((ss, qs))
        return Hits(qs[order], ss[order], len(self), len(subj))

    @flags("ignore_strand")
    def distance(self, other: GRanges, ignore_strand: bool = False) -> np.ma.MaskedArray:
        """Parallel gap; NA when seqnames differ or strands are "+" vs "-" (unless ignored)."""
        d = self._ranges.distance(other._ranges)
        n = len(d)
        if n == 0:
            return d
        reps = np.arange(n)
        ia, ib = reps % len(self), reps % len(other)
        qa, sb = self._union_codes(other)
        bad = qa[ia] != sb[ib]
        if not ignore_strand:
            bad |= (self._strand[ia].astype(int) + other._strand[ib].astype(int)) == 1  # "+" with "-"
        return np.ma.masked_array(np.ma.getdata(d), mask=np.ma.getmaskarray(d) | bad)

    # -- coverage -------------------------------------------------------------------
    def coverage(self, shift: Any = 0, width: Any = None, weight: Any = 1) -> dict[str, Rle]:
        """Per-seqlevel coverage (strand ignored); width defaults to seqlengths, else the max end."""
        if isinstance(weight, str):
            if self._mcols is None or weight not in self._mcols.columns:
                raise ValueError(f"no mcols column named {weight!r}")
            weight = self._mcols[weight].to_numpy()
        start, wt = self._ranges._coverage_inputs(shift, weight)
        levels = self._seqinfo.seqnames
        widths = self._per_seqlevel(width, "width") if width is not None else list(self._seqinfo.seqlengths)
        if np.ma.isMaskedArray(width) or (width is not None and not isinstance(width, (int, np.integer, dict, list, tuple, np.ndarray))):
            raise ValueError("'width' must be NULL or an integer vector")
        out: dict[str, Rle] = {}
        end = start + self.width - 1
        for code, name in enumerate(levels):
            sel = self._seq == code
            total = widths[code]
            if total is None:
                total = int(max(end[sel].max(initial=0), 0))
            if total < 0:
                raise ValueError("'width' cannot be negative")
            st, wd, w_ = start[sel], self.width[sel], wt[sel]
            circ_len = self._seqinfo.seqlengths[code] if self._seqinfo.is_circular[code] else None
            if circ_len:
                if width is not None and total > circ_len:
                    raise ValueError(f"'width' cannot be greater than the length of the underlying circular sequence")
                st, wd, w_ = self._wrap_circular(st, wd, w_, circ_len)
            vals, lens = coverage_kernel.coverage(st, wd, w_, int(total))
            out[name] = Rle(vals, lens)
        return out

    @staticmethod
    def _wrap_circular(
        start: np.ndarray, width: np.ndarray, weight: np.ndarray, length: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Fold ranges onto a circular sequence of `length`: whole turns, then up to two pieces."""
        ok = width > 0
        start, width, weight = start[ok], width[ok], weight[ok]
        first = (start - 1) % length + 1
        turns, rem = width // length, width % length
        spill = np.maximum(first + rem - 1 - length, 0)
        st = np.concatenate([np.ones(len(start), dtype=np.int64), first, np.ones(len(start), dtype=np.int64)])
        wd = np.concatenate([np.where(turns > 0, length, 0), np.minimum(rem, length - first + 1), spill])
        wt = np.concatenate([weight * turns, weight, weight])
        return st, wd, wt

    # -- ordering --------------------------------------------------------------------
    def _order_keys(self, ignore_strand: bool = False) -> list[np.ndarray]:
        keys = [self._seq.astype(np.int64)]
        if not ignore_strand:
            keys.append(self._strand.astype(np.int64))
        return [*keys, self.start, self.width]

    @flags("decreasing")
    def order(self, decreasing: bool = False) -> np.ndarray:
        return order_kernel.order(self._order_keys(), decreasing)

    @flags("decreasing", "ignore_strand")
    def sort(self, decreasing: bool = False, ignore_strand: bool = False) -> Self:
        return self._take(order_kernel.order(self._order_keys(ignore_strand), decreasing))

    @flags("ignore_strand")
    def rank(self, ties_method: str = "average", ignore_strand: bool = False) -> np.ndarray:
        return order_kernel.rank(self._order_keys(ignore_strand), _args.match_arg(ties_method, order_kernel.TIES, "ties.method"))

    @flags("from_last")
    def duplicated(self, from_last: bool = False) -> np.ndarray:
        return order_kernel.duplicated(self._order_keys(), from_last)

    @flags("from_last")
    def unique(self, from_last: bool = False) -> Self:
        return self[~self.duplicated(from_last)]
