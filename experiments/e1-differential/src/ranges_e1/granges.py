"""GRanges: IRanges plus seqnames, strand and seqinfo."""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _meta
from ._kernels import inter, intra, setops
from ._meta import Names
from .iranges import IRanges
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
        missing = [s for s in sn if s not in index]
        if missing:
            raise ValueError("'seqnames' contains sequence names with no entries in 'seqinfo'")
        codes = np.array([index[s] for s in sn], dtype=np.int32)
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

    def shift(self, shift: Any = 0, use_names: bool = True) -> Self:
        """Shift along the genome; strand is ignored (R)."""
        return self._moved(self._ranges.shift(shift), use_names)

    def narrow(self, start: Any = None, end: Any = None, width: Any = None, use_names: bool = True) -> Self:
        """Narrow relative to each range's left end; strand is ignored (R)."""
        return self._moved(self._ranges.narrow(start, end, width), use_names)

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

    def range(self, *others: GRanges, ignore_strand: bool = False) -> Self:
        """Smallest range per (seqname, strand) group present."""
        x = self.concat([self, *others]) if others else self
        red = inter.range_(x._gid(ignore_strand), x.start, x.width)
        return x._from_groups(red.group, red.start, red.width)

    def reduce(
        self, drop_empty_ranges: bool = False, min_gapwidth: int = 1, with_revmap: bool = False,
        ignore_strand: bool = False,
    ) -> Self:
        if min_gapwidth < 0:
            raise ValueError("'min.gapwidth' must be non-negative")
        red = inter.reduce(self._gid(ignore_strand), self.start, self.width, int(min_gapwidth), bool(drop_empty_ranges))
        mc = _meta.revmap_frame(red.members) if with_revmap else None
        return self._from_groups(red.group, red.start, red.width, mc)

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
        return [None if (x is None or (isinstance(x, float) and x != x)) else int(x) for x in vals]

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

    def union(self, other: GRanges, ignore_strand: bool = False) -> Self:
        return self._set_op(other, setops.union, ignore_strand)

    def intersect(self, other: GRanges, ignore_strand: bool = False) -> Self:
        return self._set_op(other, setops.intersect, ignore_strand)

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

    def _merged_strand(self, st_b: np.ndarray, ok: np.ndarray, ignore_strand: bool) -> np.ndarray:
        """Strand of a compatible pair: x's, unless x is "*" (then y's); x's when strands are ignored."""
        if ignore_strand:
            return self._strand
        return np.where(ok & (self._strand == STAR), st_b, self._strand).astype(np.int8)

    def punion(self, other: GRanges, fill_gap: bool = False, ignore_strand: bool = False) -> Self:
        si, ys, yw, st_b, ok = self._pair_with(other, recycle_other=False, ignore_strand=ignore_strand)
        if not ok.all():
            raise ValueError("'x' and 'y' elements must have compatible 'seqnames' and 'strand' values")
        s, w = setops.punion(self.start, self.width, ys, yw, bool(fill_gap))
        return self._new(self._seq_in(si), s, w, self._strand, self._names, None, si, check=True)

    def _seq_in(self, si: Seqinfo) -> np.ndarray:
        """This object's seqlevel codes expressed in the (merged) seqinfo `si`."""
        index = {s: i for i, s in enumerate(si.seqnames)}
        return np.array([index[s] for s in self._seqinfo.seqnames], dtype=np.int32)[self._seq]

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
        strand = np.where(hit, self._merged_strand(st_b, ok, ignore_strand), self._strand).astype(np.int8)
        names, mcols = self._names, self._mcols
        if drop_nohit_ranges:
            keep = np.flatnonzero(hit)
            return self._new(
                self._seq_in(si)[keep], start[keep], width[keep], strand[keep],
                _meta.take_names(names, keep), _meta.take_mcols(mcols, keep), si, check=True,
            )
        hit_col = pl.DataFrame([pl.Series("hit", hit.tolist(), dtype=pl.Boolean)])
        mc = hit_col if mcols is None else pl.concat([mcols, hit_col], how="horizontal")
        return self._new(self._seq_in(si), start, width, strand, names, mc, si, check=True)

    def psetdiff(self, other: GRanges, ignore_strand: bool = False) -> Self:
        si, ys, yw, st_b, ok = self._pair_with(other, recycle_other=False, ignore_strand=ignore_strand)
        start, width = self.start.copy(), self.width.copy()
        idx = np.flatnonzero(ok)
        if len(idx):
            s, w = setops.psetdiff(start[idx], width[idx], ys[idx], yw[idx])
            start[idx], width[idx] = s, w
        strand = np.where(ok, self._merged_strand(st_b, ok, ignore_strand), self._strand).astype(np.int8)
        return self._new(self._seq_in(si), start, width, strand, self._names, None, si, check=True)
