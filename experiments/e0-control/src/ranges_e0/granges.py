"""GRanges: IRanges + seqnames + strand + seqinfo.

Design: a ``GRanges`` stores ``start``/``width`` int64 arrays plus integer codes
for the seqname (index into ``seqinfo``) and the strand (0 ``+``, 1 ``-``,
2 ``*``, which is R's level order, so sorting by code sorts like R). Every
algorithm is the IRanges kernel applied per (seqlevel, strand) group.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import polars as pl

from . import _coverage, _inter, _intra, _nearest, _order, _overlaps, _setops, _validate as V
from .hits import Hits
from .iranges import NO_MATCH, IRanges
from .rle import Rle
from .seqinfo import Seqinfo

I64 = np.int64
STRAND_LABELS = np.array(["+", "-", "*"], dtype=object)
_STRAND_CODE = {"+": 0, "-": 1, "*": 2}
PLUS, MINUS, STAR = 0, 1, 2


class OutOfBoundWarning(UserWarning):
    """Ranges extend past the known length of their (non-circular) sequence."""


class GRanges:
    """Genomic ranges. Immutable by convention; operations return new objects."""

    __slots__ = ("_seq", "_start", "_width", "_strand", "_seqinfo", "_names", "_mcols")

    def __init__(
        self,
        seqnames: str | Sequence[str] | np.ndarray,
        ranges: IRanges,
        strand: str | Sequence[str] | None = None,
        mcols: pl.DataFrame | dict[str, Sequence[Any]] | None = None,
        seqinfo: Seqinfo | None = None,
        seqlengths: Mapping[str, int | None] | None = None,
        names: Sequence[str] | None = None,
    ) -> None:
        sn = [seqnames] if isinstance(seqnames, str) else [str(s) for s in seqnames]
        n = max(len(sn), len(ranges))
        if len(sn) == 0 or len(ranges) == 0:
            n = 0
        if len(sn) != n and len(sn) != 1:
            raise ValueError(
                f"'seqnames' must have the length of the object to construct ({n}) or length 1"
            )
        if len(ranges) != n and len(ranges) != 1:
            raise ValueError(f"'ranges' must have the length of the object to construct ({n}) or length 1")
        if len(sn) == 1 and n != 1:
            sn = sn * n
        if len(ranges) == 1 and n != 1:
            ranges = IRanges._make(
                np.repeat(ranges._start, n), np.repeat(ranges._width, n), None
            )
        if seqinfo is not None and seqlengths is not None:
            raise ValueError("supply only one of 'seqinfo' and 'seqlengths'")
        if seqinfo is None:
            if seqlengths is not None:
                seqinfo = Seqinfo.from_mapping(seqlengths)
            else:
                seqinfo = Seqinfo(list(dict.fromkeys(sn)))
        missing = sorted({s for s in sn if s not in seqinfo})
        if missing:
            raise ValueError("'seqnames' contains sequence names with no entries in 'seqinfo': " + ", ".join(missing))
        self._seq = np.array([seqinfo.index(s) for s in sn], dtype=np.int32)
        self._strand = _as_strand_codes(strand, n)
        self._start = ranges._start
        self._width = ranges._width
        self._seqinfo = seqinfo
        self._names = V.as_names(names, n) if names is not None else ranges._names
        self._mcols = V.as_mcols(mcols, n)
        self._warn_out_of_bound()

    @classmethod
    def _make(
        cls,
        seq: np.ndarray,
        start: np.ndarray,
        width: np.ndarray,
        strand: np.ndarray,
        seqinfo: Seqinfo,
        names: np.ndarray | None = None,
        mcols: pl.DataFrame | None = None,
        check_bounds: bool = False,
    ) -> GRanges:
        """Internal constructor from validated arrays."""
        if len(width) and width.min() < 0:
            raise ValueError("negative width")
        obj = object.__new__(cls)
        obj._seq = seq.astype(np.int32, copy=False)
        obj._start = start
        obj._width = width
        obj._strand = strand.astype(np.int8, copy=False)
        obj._seqinfo = seqinfo
        obj._names = names
        obj._mcols = pl.DataFrame() if mcols is None else mcols
        if check_bounds:
            obj._warn_out_of_bound()
        return obj

    @classmethod
    def empty(cls, seqinfo: Seqinfo | None = None) -> GRanges:
        z = np.array([], dtype=I64)
        return cls._make(z.astype(np.int32), z, z, z.astype(np.int8), seqinfo or Seqinfo([]))

    # ----- accessors -------------------------------------------------------
    def __len__(self) -> int:
        return len(self._start)

    @property
    def seqnames(self) -> np.ndarray:
        """Sequence name of each range (object array of ``str``)."""
        levels = np.array(self._seqinfo.seqnames, dtype=object)
        return levels[self._seq] if len(self._seq) else np.array([], dtype=object)

    @property
    def strand(self) -> np.ndarray:
        """Strand of each range: ``"+"``, ``"-"`` or ``"*"`` (object array)."""
        return STRAND_LABELS[self._strand]

    @property
    def start(self) -> np.ndarray:
        return _ro(self._start)

    @property
    def width(self) -> np.ndarray:
        return _ro(self._width)

    @property
    def end(self) -> np.ndarray:
        return self._start + self._width - 1

    @property
    def ranges(self) -> IRanges:
        """The ranges as an ``IRanges`` (carrying this object's names)."""
        return IRanges._make(self._start, self._width, self._names)

    @property
    def names(self) -> list[str] | None:
        return None if self._names is None else list(self._names)

    @property
    def mcols(self) -> pl.DataFrame:
        return self._mcols

    @property
    def seqinfo(self) -> Seqinfo:
        return self._seqinfo

    @property
    def seqlevels(self) -> tuple[str, ...]:
        return self._seqinfo.seqnames

    @property
    def seqlengths(self) -> dict[str, int | None]:
        return self._seqinfo.seqlengths

    @property
    def is_circular(self) -> dict[str, bool | None]:
        return self._seqinfo.is_circular

    @property
    def genome(self) -> dict[str, str | None]:
        return self._seqinfo.genome

    # ----- copy-with helpers -----------------------------------------------
    def with_names(self, names: Sequence[str] | None) -> GRanges:
        return self._replace(names=V.as_names(names, len(self)), set_names=True)

    def with_mcols(self, mcols: pl.DataFrame | dict[str, Sequence[Any]] | None) -> GRanges:
        return self._replace(mcols=V.as_mcols(mcols, len(self)))

    def with_strand(self, strand: str | Sequence[str] | None) -> GRanges:
        """R ``strand<-``."""
        return self._replace(strand=_as_strand_codes(strand, len(self)))

    def with_seqinfo(self, seqinfo: Seqinfo) -> GRanges:
        """Replace the seqinfo, re-coding seqnames by name. Every sequence used
        by a range must exist in the new seqinfo (R ``seqinfo<-`` is stricter
        about dropping levels; here that is an error)."""
        used = np.unique(self._seq)
        old = self._seqinfo.seqnames
        lookup = np.full(len(old), -1, dtype=np.int32)
        for i in used:
            if old[i] not in seqinfo:
                raise ValueError(f"sequence {old[i]!r} is used by ranges but absent from the new seqinfo")
        for i, name in enumerate(old):
            if name in seqinfo:
                lookup[i] = seqinfo.index(name)
        return GRanges._make(
            lookup[self._seq] if len(self) else self._seq,
            self._start, self._width, self._strand, seqinfo, self._names, self._mcols, check_bounds=True,
        )

    def _replace(
        self,
        *,
        seq: np.ndarray | None = None,
        start: np.ndarray | None = None,
        width: np.ndarray | None = None,
        strand: np.ndarray | None = None,
        seqinfo: Seqinfo | None = None,
        names: np.ndarray | None = None,
        set_names: bool = False,
        mcols: pl.DataFrame | None = None,
        check_bounds: bool = False,
    ) -> GRanges:
        return GRanges._make(
            self._seq if seq is None else seq,
            self._start if start is None else start,
            self._width if width is None else width,
            self._strand if strand is None else strand,
            self._seqinfo if seqinfo is None else seqinfo,
            names if set_names else self._names,
            self._mcols if mcols is None else mcols,
            check_bounds,
        )

    # ----- bounds ------------------------------------------------------------
    def _out_of_bound(self) -> np.ndarray:
        """Bool per range: extends past a known, non-circular sequence length."""
        si = self._seqinfo
        lens = np.array([-1 if l is None or c else l for l, c in zip(si._lengths, si._circular)], dtype=I64)
        if len(self) == 0:
            return np.zeros(0, dtype=bool)
        L = lens[self._seq]
        known = L >= 0
        return known & ((self._start < 1) | (self._start + self._width - 1 > L))

    def _warn_out_of_bound(self) -> None:
        oob = self._out_of_bound()
        if oob.any():
            levels = sorted({self._seqinfo.seqnames[i] for i in self._seq[oob]})
            warnings.warn(
                f"GRanges object contains {int(oob.sum())} out-of-bound range(s) located on "
                f"sequence(s) {', '.join(levels)}; use trim() to trim them",
                OutOfBoundWarning,
                stacklevel=3,
            )

    def trim(self) -> GRanges:
        """R ``trim``: clip ranges to ``[1, seqlength]`` on non-circular
        sequences of known length. A range entirely beyond the end becomes the
        empty range ``(len + 1, len)``; entirely before 1 becomes ``(1, 0)``."""
        si = self._seqinfo
        lens = np.array([-1 if l is None or c else l for l, c in zip(si._lengths, si._circular)], dtype=I64)
        if len(self) == 0:
            return self
        L = lens[self._seq]
        known = L >= 0
        end = self._start + self._width - 1
        new_start = np.where(known, np.clip(self._start, 1, np.where(known, L, 0) + 1), self._start)
        new_end = np.where(known, np.clip(end, 0, np.where(known, L, 0)), end)
        return self._replace(start=new_start, width=new_end - new_start + 1)

    # ----- subsetting / combining --------------------------------------------
    def __getitem__(self, index: Any) -> GRanges:
        return self._take(V.resolve_index(index, len(self), self._names))

    def _take(self, idx: np.ndarray) -> GRanges:
        return GRanges._make(
            self._seq[idx], self._start[idx], self._width[idx], self._strand[idx], self._seqinfo,
            None if self._names is None else self._names[idx], V.take_rows(self._mcols, idx),
        )

    @classmethod
    def concat(cls, *objs: GRanges) -> GRanges:
        """R ``c()``: merges seqinfo (conflicts raise ``ValueError``), pads names
        with ``""`` and aligns mcols by column name."""
        if not objs or not all(isinstance(o, GRanges) for o in objs):
            raise TypeError("can only concatenate GRanges with GRanges")
        si = objs[0]._seqinfo
        for o in objs[1:]:
            si = si.merge(o._seqinfo)
        seq = np.concatenate([_recode(o._seq, o._seqinfo, si) for o in objs])
        names = None
        if any(o._names is not None for o in objs):
            names = np.concatenate(
                [o._names if o._names is not None else np.full(len(o), "", dtype=object) for o in objs]
            )
        return GRanges._make(
            seq,
            np.concatenate([o._start for o in objs]),
            np.concatenate([o._width for o in objs]),
            np.concatenate([o._strand for o in objs]),
            si,
            names,
            V.concat_mcols([o._mcols for o in objs], [len(o) for o in objs]),
            check_bounds=True,
        )

    # ----- intra-range ------------------------------------------------------
    def _neg(self, ignore_strand: bool) -> np.ndarray:
        """Rows treated as minus strand (``*`` behaves like ``+``)."""
        return np.zeros(len(self), dtype=bool) if ignore_strand else self._strand == MINUS

    def _intra_result(self, start: np.ndarray, width: np.ndarray, use_names: bool) -> GRanges:
        return self._replace(
            start=start, width=width, names=None if use_names else None, set_names=not use_names,
            check_bounds=True,
        )

    def shift(self, shift: Any, use_names: bool = True) -> GRanges:
        """R ``shift``: not strand-aware."""
        return self._intra_result(*_intra.shift_k(self._start, self._width, shift), use_names)

    def narrow(self, start: Any = None, end: Any = None, width: Any = None, use_names: bool = True) -> GRanges:
        """R ``narrow``: offsets relative to each range, not strand-aware."""
        return self._intra_result(*_intra.narrow_k(self._start, self._width, start, end, width), use_names)

    def resize(
        self,
        width: Any,
        fix: str | Sequence[str] = "start",
        use_names: bool = True,
        ignore_strand: bool = False,
    ) -> GRanges:
        """R ``resize``: ``fix="start"`` anchors the *transcription* start, i.e.
        the end of minus-strand ranges (``*`` behaves like ``+``)."""
        fixa = _intra.as_fix(fix, len(self))
        neg = self._neg(ignore_strand)
        flipped = np.where(fixa == "start", "end", np.where(fixa == "end", "start", fixa))
        fixa = np.where(neg, flipped, fixa)
        return self._intra_result(*_intra.resize_k(self._start, self._width, width, fixa), use_names)

    def flank(
        self,
        width: Any,
        start: bool | Sequence[bool] = True,
        both: bool = False,
        use_names: bool = True,
        ignore_strand: bool = False,
    ) -> GRanges:
        """R ``flank``: ``start=True`` is the upstream side (the end of
        minus-strand ranges)."""
        at = np.asarray(start, dtype=bool)
        at = np.repeat(at, len(self)) if at.ndim == 0 or len(at) == 1 else at
        at = np.where(self._neg(ignore_strand), ~at, at)
        return self._intra_result(*_intra.flank_k(self._start, self._width, width, at, both), use_names)

    def promoters(self, upstream: Any = 2000, downstream: Any = 200, use_names: bool = True) -> GRanges:
        """R ``promoters``: window around the TSS (start of ``+``/``*``, end of ``-`` ranges)."""
        s, w = _intra.promoters_k(self._start, self._width, upstream, downstream, self._neg(False))
        return self._intra_result(s, w, use_names)

    # ----- grouping seam ----------------------------------------------------
    def _group_key(self, ignore_strand: bool) -> np.ndarray:
        strand = np.full(len(self), STAR, dtype=I64) if ignore_strand else self._strand.astype(I64)
        return self._seq.astype(I64) * 3 + strand

    def _groups(self, ignore_strand: bool) -> dict[int, np.ndarray]:
        """``{seq * 3 + strand: row indices}`` in key order (seqlevel, then ``+ - *``).
        With ``ignore_strand`` every row is in the ``*`` group of its seqlevel."""
        key = self._group_key(ignore_strand)
        if len(key) == 0:
            return {}
        order = np.argsort(key, kind="stable")
        ks = key[order]
        cuts = np.flatnonzero(np.diff(ks)) + 1
        return {int(ks[i]): part for i, part in zip(np.concatenate(([0], cuts)), np.split(order, cuts))}

    def _assemble(
        self,
        parts: list[tuple[int, np.ndarray, np.ndarray]],
        seqinfo: Seqinfo | None = None,
        revmap: list[np.ndarray] | None = None,
    ) -> GRanges:
        """Build a GRanges from per-group results ``(key, start, width)``."""
        seqinfo = self._seqinfo if seqinfo is None else seqinfo
        if parts:
            seq = np.concatenate([np.full(len(st), k // 3, dtype=np.int32) for k, st, _ in parts])
            strand = np.concatenate([np.full(len(st), k % 3, dtype=np.int8) for k, st, _ in parts])
            start = np.concatenate([st for _, st, _ in parts])
            width = np.concatenate([w for _, _, w in parts])
        else:
            seq, strand = np.array([], dtype=np.int32), np.array([], dtype=np.int8)
            start = width = np.array([], dtype=I64)
        mcols = None
        if revmap is not None:
            mcols = pl.DataFrame({"revmap": pl.Series("revmap", [r.tolist() for r in revmap], dtype=pl.List(pl.Int64))})
        return GRanges._make(seq, start, width, strand, seqinfo, None, mcols)

    # ----- inter-range -------------------------------------------------------
    def range(self, *others: GRanges, ignore_strand: bool = False) -> GRanges:
        """R ``range``: one range per (seqlevel, strand) group present."""
        pooled = GRanges.concat(self, *others) if others else self
        parts = []
        for key, idx in pooled._groups(ignore_strand).items():
            st, w = _inter.range_k(pooled._start[idx], pooled._width[idx])
            parts.append((key, st, w))
        return pooled._assemble(parts)

    def reduce(
        self,
        min_gapwidth: int = 1,
        drop_empty_ranges: bool = False,
        with_revmap: bool = False,
        ignore_strand: bool = False,
    ) -> GRanges:
        """R ``reduce``, per (seqlevel, strand) group; ``revmap`` holds 0-based
        indices into this object."""
        parts: list[tuple[int, np.ndarray, np.ndarray]] = []
        revs: list[np.ndarray] = []
        for key, idx in self._groups(ignore_strand).items():
            st, w, rev = _inter.reduce_k(
                self._start[idx], self._width[idx], min_gapwidth, drop_empty_ranges, with_revmap
            )
            parts.append((key, st, w))
            if rev is not None:
                revs.extend(idx[r] for r in rev)
        return self._assemble(parts, revmap=revs if with_revmap else None)

    def disjoin(self, with_revmap: bool = False, ignore_strand: bool = False) -> GRanges:
        """R ``disjoin``, per (seqlevel, strand) group."""
        parts: list[tuple[int, np.ndarray, np.ndarray]] = []
        revs: list[np.ndarray] = []
        for key, idx in self._groups(ignore_strand).items():
            st, w, rev = _inter.disjoin_k(self._start[idx], self._width[idx], with_revmap)
            parts.append((key, st, w))
            if rev is not None:
                revs.extend(idx[r] for r in rev)
        return self._assemble(parts, revmap=revs if with_revmap else None)

    def gaps(self, start: int | None = 1, end: int | None = None, ignore_strand: bool = False) -> GRanges:
        """R ``gaps``: uncovered positions in ``[start, end]`` for *every* seqlevel
        and strand (``+ - *``, or just ``*`` with ``ignore_strand``), including
        groups without ranges. ``end=None`` uses the seqlengths; a sequence of
        unknown length has no trailing gap and no gaps at all when empty."""
        groups = self._groups(ignore_strand)
        strands = (STAR,) if ignore_strand else (PLUS, MINUS, STAR)
        empty = np.array([], dtype=I64)
        parts = []
        for i, length in enumerate(self._seqinfo._lengths):
            gap_end = end if end is not None else length
            for strand in strands:
                key = i * 3 + strand
                if ignore_strand:
                    idx = np.concatenate([v for k, v in groups.items() if k // 3 == i] or [empty])
                else:
                    idx = groups.get(key, empty)
                st, w = _inter.gaps_k(self._start[idx], self._width[idx], start, gap_end)
                if len(st):
                    parts.append((key, st, w))
        return self._assemble(parts)

    # ----- set operations -----------------------------------------------------
    def _pair_groups(
        self, other: GRanges, ignore_strand: bool
    ) -> tuple[Seqinfo, GRanges, GRanges]:
        """Both operands with seqinfo merged and seq codes re-expressed in it."""
        si = self._seqinfo.merge(other._seqinfo)
        x = self._replace(seq=_recode(self._seq, self._seqinfo, si), seqinfo=si)
        y = other._replace(seq=_recode(other._seq, other._seqinfo, si), seqinfo=si)
        return si, x, y

    def _setop(self, other: GRanges, kernel: Any, keys_from: str, ignore_strand: bool) -> GRanges:
        si, x, y = self._pair_groups(other, ignore_strand)
        gx, gy = x._groups(ignore_strand), y._groups(ignore_strand)
        empty = np.array([], dtype=I64)
        if keys_from == "any":
            keys = sorted(set(gx) | set(gy))
        elif keys_from == "x":
            keys = sorted(gx)
        else:
            keys = sorted(set(gx) & set(gy))
        parts = []
        for key in keys:
            ix, iy = gx.get(key, empty), gy.get(key, empty)
            st, w = kernel(x._start[ix], x._width[ix], y._start[iy], y._width[iy])
            if len(st):
                parts.append((key, st, w))
        return x._assemble(parts, si)

    def union(self, other: GRanges, ignore_strand: bool = False) -> GRanges:
        """R ``union``: per (seqlevel, strand) group; ``*`` is its own group."""
        return self._setop(other, _setops.union_k, "any", ignore_strand)

    def intersect(self, other: GRanges, ignore_strand: bool = False) -> GRanges:
        return self._setop(other, _setops.intersect_k, "common", ignore_strand)

    def setdiff(self, other: GRanges, ignore_strand: bool = False) -> GRanges:
        return self._setop(other, _setops.setdiff_k, "x", ignore_strand)

    def _pairwise_setup(self, other: GRanges, ignore_strand: bool) -> tuple[Seqinfo, GRanges, GRanges, np.ndarray]:
        if len(self) != len(other):
            raise ValueError("'x' and 'y' must have the same length")
        si, x, y = self._pair_groups(other, ignore_strand)
        if ignore_strand:
            y = y._replace(strand=x._strand)
        compat = (x._seq == y._seq) & ((x._strand == STAR) | (y._strand == STAR) | (x._strand == y._strand))
        return si, x, y, compat

    def punion(self, other: GRanges, fill_gap: bool = False, ignore_strand: bool = False) -> GRanges:
        """R ``punion``: strand/seqname must be compatible (``*`` matches anything);
        the result keeps this object's strand and names, drops mcols."""
        si, x, y, compat = self._pairwise_setup(other, ignore_strand)
        if not compat.all():
            raise ValueError("'x' and 'y' elements must have compatible 'seqnames' and 'strand' values")
        s, w = _setops.punion_k(x._start, x._width, y._start, y._width, fill_gap)
        return x._replace(start=s, width=w, mcols=pl.DataFrame())

    def pintersect(
        self,
        other: GRanges,
        drop_nohit_ranges: bool = False,
        ignore_strand: bool = False,
        strict_strand: bool = False,
    ) -> GRanges:
        """R ``pintersect``. Pairs on different seqnames/strands (or disjoint) are
        "no hit": they become the empty range at ``start(x)`` and the result gets a
        logical mcol ``hit`` (unless ``drop_nohit_ranges``, which removes them).
        Unless ``strict_strand``, ``*`` on one side adopts the other side's strand."""
        if len(other) != len(self):
            if len(other) != 1:
                raise ValueError("'y' must have the length of 'x' or length 1")
            other = other._take(np.zeros(len(self), dtype=np.int64))
        x_orig_strand = self._strand
        si, x, y = self._pair_groups(other, ignore_strand)
        if ignore_strand:
            y = y._replace(strand=x._strand)
        elif not strict_strand:
            xs, ys = x._strand.copy(), y._strand.copy()
            xs = np.where((x._strand == STAR) & (y._strand != STAR), y._strand, xs)
            ys = np.where((y._strand == STAR) & (x._strand != STAR), x._strand, ys)
            x, y = x._replace(strand=xs), y._replace(strand=ys)
        new_start = np.maximum(x._start, y._start)
        new_end = np.minimum(x.end, y.end)
        is_hit = (x._seq == y._seq) & (x._strand == y._strand) & (new_end >= new_start - 1)
        if drop_nohit_ranges:
            keep = np.flatnonzero(is_hit)
            x = x._take(keep)
            return x._replace(start=new_start[keep], width=new_end[keep] - new_start[keep] + 1)
        nohit = ~is_hit
        new_start = np.where(nohit, x._start, new_start)
        new_end = np.where(nohit, x._start - 1, new_end)
        strand = x._strand
        if not (ignore_strand or strict_strand):
            strand = np.where(nohit, x_orig_strand, strand)
        mcols = x._mcols.with_columns(pl.Series("hit", is_hit.tolist(), dtype=pl.Boolean)) if x._mcols.width else (
            pl.DataFrame({"hit": pl.Series("hit", is_hit.tolist(), dtype=pl.Boolean)})
        )
        return x._replace(start=new_start, width=new_end - new_start + 1, strand=strand, mcols=mcols)

    def psetdiff(self, other: GRanges, ignore_strand: bool = False) -> GRanges:
        """R ``psetdiff``: pairs that are not compatible keep ``x``'s range; a ``*``
        strand in ``x`` is resolved from ``y``. Drops mcols."""
        si, x, y, ok = self._pairwise_setup(other, ignore_strand)
        s, w = x._start.copy(), x._width.copy()
        if ok.any():
            ps, pw = _setops.psetdiff_k(x._start[ok], x._width[ok], y._start[ok], y._width[ok])
            s[ok], w[ok] = ps, pw
        strand = np.where(x._strand == STAR, y._strand, x._strand)
        return x._replace(start=s, width=w, strand=strand, mcols=pl.DataFrame())

    # ----- overlaps -----------------------------------------------------------
    def _align_seq(self, subject: GRanges) -> np.ndarray:
        """Subject seqname codes expressed in *this* object's seqinfo (-1: level absent here)."""
        lookup = np.array(
            [self._seqinfo.index(n) if n in self._seqinfo else -1 for n in subject._seqinfo.seqnames],
            dtype=np.int64,
        )
        return lookup[subject._seq] if len(subject) else np.array([], dtype=np.int64)

    def _overlap_pairs(
        self,
        subject: GRanges,
        maxgap: int,
        minoverlap: int,
        type: str,
        ignore_strand: bool,
    ) -> tuple[np.ndarray, np.ndarray]:
        """(query, subject) index pairs on the same seqname with compatible strands
        (``*`` matches any strand), sorted by query then subject."""
        _overlaps.check_args(type, maxgap, minoverlap)
        s_seq = self._align_seq(subject)
        qi_parts: list[np.ndarray] = []
        si_parts: list[np.ndarray] = []
        for seq in np.unique(self._seq):
            q_rows = np.flatnonzero(self._seq == seq)
            s_rows = np.flatnonzero(s_seq == seq)
            if len(s_rows) == 0:
                continue
            for strand in ((STAR,) if ignore_strand else (PLUS, MINUS, STAR)):
                if ignore_strand:
                    qr, sr = q_rows, s_rows
                else:
                    qr = q_rows[self._strand[q_rows] == strand]
                    sr = s_rows if strand == STAR else s_rows[subject._strand[s_rows] != (1 - strand)]
                if len(qr) == 0 or len(sr) == 0:
                    continue
                qi, si = _overlaps.overlap_pairs(
                    self._start[qr], self._width[qr], subject._start[sr], subject._width[sr],
                    type, maxgap, minoverlap,
                )
                qi_parts.append(qr[qi])
                si_parts.append(sr[si])
        if not qi_parts:
            return np.array([], dtype=I64), np.array([], dtype=I64)
        qi, si = np.concatenate(qi_parts), np.concatenate(si_parts)
        o = np.lexsort((si, qi))
        return qi[o], si[o]

    def find_overlaps(
        self,
        subject: GRanges | None = None,
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",
        select: str = "all",
        ignore_strand: bool = False,
        drop_self: bool = False,
        drop_redundant: bool = False,
    ) -> Hits | np.ndarray:
        """R ``findOverlaps`` for GRanges: same seqname, strands compatible
        (``*`` overlaps ``+`` and ``-``) unless ``ignore_strand``; the IRanges
        predicate for ``type``/``maxgap``/``minoverlap`` (see IRanges.find_overlaps)."""
        if subject is None:
            qi, si = self._overlap_pairs(self, maxgap, minoverlap, type, ignore_strand)
            keep = np.ones(len(qi), dtype=bool)
            if drop_self:
                keep &= qi != si
            if drop_redundant:
                keep &= qi <= si
            qi, si = qi[keep], si[keep]
            nsub = len(self)
        else:
            if drop_self or drop_redundant:
                raise ValueError("'drop_self' and 'drop_redundant' require 'subject' to be omitted")
            qi, si = self._overlap_pairs(subject, maxgap, minoverlap, type, ignore_strand)
            nsub = len(subject)
        if select == "all":
            return Hits(qi, si, len(self), nsub)
        return _overlaps.select_hits(qi, si, len(self), select, NO_MATCH)

    def count_overlaps(
        self, subject: GRanges, maxgap: int = -1, minoverlap: int = 0, type: str = "any",
        ignore_strand: bool = False,
    ) -> np.ndarray:
        qi, _ = self._overlap_pairs(subject, maxgap, minoverlap, type, ignore_strand)
        return np.bincount(qi, minlength=len(self)).astype(np.int64)

    def overlaps_any(
        self, subject: GRanges, maxgap: int = -1, minoverlap: int = 0, type: str = "any",
        ignore_strand: bool = False,
    ) -> np.ndarray:
        return self.count_overlaps(subject, maxgap, minoverlap, type, ignore_strand) > 0

    def subset_by_overlaps(
        self, ranges: GRanges, maxgap: int = -1, minoverlap: int = 0, type: str = "any",
        invert: bool = False, ignore_strand: bool = False,
    ) -> GRanges:
        hit = self.overlaps_any(ranges, maxgap, minoverlap, type, ignore_strand)
        return self._take(np.flatnonzero(~hit if invert else hit))

    # ----- nearest ---------------------------------------------------------------
    def _precede_follow_pairs(
        self, subject: GRanges, where: str, ignore_strand: bool
    ) -> tuple[np.ndarray, np.ndarray]:
        """All tied (query, subject) pairs for ``precede`` / ``follow``.

        A ``+`` query looks right for precede (left for follow) among ``+``/``*``
        subjects; a ``-`` query looks the other way among ``-``/``*`` subjects; a
        ``*`` query does both (``+``/``*`` subjects going the + way, only ``-``
        subjects going the - way) and keeps the closer side(s). Same seqname only."""
        right_for_plus = where == "precede"
        s_seq = self._align_seq(subject)
        xstr = np.full(len(self), PLUS) if ignore_strand else self._strand
        sstr = np.full(len(subject), PLUS) if ignore_strand else subject._strand
        qp: list[np.ndarray] = []
        sp: list[np.ndarray] = []

        def directional(look_right: bool, qr: np.ndarray, sr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            if look_right:
                q, s = _nearest.precede_k(self._start[qr], self._width[qr], subject._start[sr], "all")
            else:
                q, s = _nearest.follow_k(self._start[qr], subject._start[sr], subject._width[sr], "all")
            return qr[q], sr[s]

        for seq in np.unique(self._seq):
            q_rows = np.flatnonzero(self._seq == seq)
            s_rows = np.flatnonzero(s_seq == seq)
            if len(s_rows) == 0:
                continue
            plus_subj = s_rows[sstr[s_rows] != MINUS]
            minus_subj = s_rows[sstr[s_rows] != PLUS]
            only_minus = s_rows[sstr[s_rows] == MINUS]
            for strand in (PLUS, MINUS, STAR):
                qr = q_rows[xstr[q_rows] == strand]
                if len(qr) == 0:
                    continue
                if strand == PLUS:
                    pieces = [directional(right_for_plus, qr, plus_subj)] if len(plus_subj) else []
                elif strand == MINUS:
                    pieces = [directional(not right_for_plus, qr, minus_subj)] if len(minus_subj) else []
                else:
                    a = directional(right_for_plus, qr, plus_subj) if len(plus_subj) else _EMPTY_PAIRS
                    b = directional(not right_for_plus, qr, only_minus) if len(only_minus) else _EMPTY_PAIRS
                    pieces = [self._closer_side(subject, a, b)]
                for q, s in pieces:
                    qp.append(q)
                    sp.append(s)
        if not qp:
            return _EMPTY_PAIRS
        q, s = np.concatenate(qp), np.concatenate(sp)
        o = np.lexsort((s, q))
        return q[o], s[o]

    def _closer_side(
        self, subject: GRanges, a: tuple[np.ndarray, np.ndarray], b: tuple[np.ndarray, np.ndarray]
    ) -> tuple[np.ndarray, np.ndarray]:
        """For ``*`` queries with hits on both sides keep only the pairs at minimum distance."""
        both = np.intersect1d(a[0], b[0])
        q = np.concatenate((a[0], b[0]))
        s = np.concatenate((a[1], b[1]))
        if len(both) == 0:
            return q, s
        d = _nearest.distance_k(self._start[q], self._width[q], subject._start[s], subject._width[s])
        in_both = np.isin(q, both)
        drop = np.zeros(len(q), dtype=bool)
        for qq in both:
            sel = q == qq
            drop |= sel & (d != d[sel].min())
        keep = ~(in_both & drop)
        return q[keep], s[keep]

    def _pf_result(self, q: np.ndarray, s: np.ndarray, n_sub: int, select: str, which: str) -> Hits | np.ndarray:
        if select == "all":
            return Hits(q, s, len(self), n_sub)
        out = np.full(len(self), NO_MATCH, dtype=I64)
        if len(q):
            # pairs are sorted by (query, subject): first / last per query
            boundary = np.flatnonzero(q[1:] != q[:-1])
            pos = np.concatenate((boundary, [len(q) - 1])) if which == "last" else np.concatenate(([0], boundary + 1))
            out[q[pos]] = s[pos]
        return out

    def precede(
        self, subject: GRanges | None = None, select: str = "first", ignore_strand: bool = False
    ) -> Hits | np.ndarray:
        """R ``precede`` for GRanges (strand-aware: on ``-`` the downstream direction
        is toward lower coordinates). Ties: lowest subject index, or all."""
        if select not in ("first", "all"):
            raise ValueError("'select' must be one of ('first', 'all')")
        sub = self if subject is None else subject
        q, s = self._precede_follow_pairs(sub, "precede", ignore_strand)
        return self._pf_result(q, s, len(sub), select, "first")

    def follow(
        self, subject: GRanges | None = None, select: str = "last", ignore_strand: bool = False
    ) -> Hits | np.ndarray:
        """R ``follow`` for GRanges. Ties: highest subject index, or all."""
        if select not in ("last", "all"):
            raise ValueError("'select' must be one of ('last', 'all')")
        sub = self if subject is None else subject
        q, s = self._precede_follow_pairs(sub, "follow", ignore_strand)
        return self._pf_result(q, s, len(sub), select, "last")

    def nearest(
        self, subject: GRanges | None = None, select: str = "arbitrary", ignore_strand: bool = False
    ) -> Hits | np.ndarray:
        """R ``nearest`` for GRanges: an overlapping/adjacent (strand-compatible)
        subject if any, else the closer of precede/follow. Equidistant: the lower
        subject index (``arbitrary``) or the precede side only (``all``), as R does."""
        if select not in ("arbitrary", "all"):
            raise ValueError("'select' must be one of ('arbitrary', 'all')")
        sub = self if subject is None else subject
        n = len(self)
        if len(sub) == 0:
            return Hits(_EMPTY_PAIRS[0], _EMPTY_PAIRS[1], n, 0) if select == "all" else np.full(n, NO_MATCH, dtype=I64)
        oq, os_ = self._overlap_pairs(sub, 0, 0, "any", ignore_strand)
        if subject is None:
            keep = oq != os_
            oq, os_ = oq[keep], os_[keep]
        has_hit = np.zeros(n, dtype=bool)
        has_hit[oq] = True
        rest = np.flatnonzero(~has_hit)
        xr = self._take(rest)
        pq, ps = xr._precede_follow_pairs(sub, "precede", ignore_strand)
        fq, fs = xr._precede_follow_pairs(sub, "follow", ignore_strand)
        p = xr._pf_result(pq, ps, len(sub), "first", "first")
        f = xr._pf_result(fq, fs, len(sub), "last", "last")
        assert isinstance(p, np.ndarray) and isinstance(f, np.ndarray)
        has_p, has_f = p >= 0, f >= 0
        pd = _nearest.distance_k(xr._start, xr._width, sub._start[np.where(has_p, p, 0)], sub._width[np.where(has_p, p, 0)]) if len(rest) else np.array([], dtype=I64)
        fd = _nearest.distance_k(xr._start, xr._width, sub._start[np.where(has_f, f, 0)], sub._width[np.where(has_f, f, 0)]) if len(rest) else np.array([], dtype=I64)
        both = has_p & has_f
        tie = both & (pd == fd)
        if select == "arbitrary":
            pnearest = np.where(both, np.where(pd == fd, p < f, pd < fd), ~has_f)
            out = _overlaps.select_hits(oq, os_, n, "first", NO_MATCH)
            out[rest] = np.where(pnearest, p, f)
            return out
        pnearest = np.where(both, pd < fd, ~has_f)
        pnearest = np.where(tie, True, pnearest)
        use_f = ~pnearest[fq] if len(fq) else np.zeros(0, dtype=bool)
        use_p = pnearest[pq] if len(pq) else np.zeros(0, dtype=bool)
        q_all = np.concatenate((oq, rest[pq[use_p]], rest[fq[use_f]]))
        s_all = np.concatenate((os_, ps[use_p], fs[use_f]))
        o = np.lexsort((s_all, q_all))
        return Hits(q_all[o], s_all[o], n, len(sub))

    def distance(self, other: GRanges, ignore_strand: bool = False) -> np.ndarray:
        """R ``distance``: gap width between paired ranges; ``NO_MATCH`` (-1) where
        the seqnames differ or (unless ``ignore_strand``) one range is ``+`` and
        the other ``-``. Length-1 arguments are recycled."""
        if len(self) == 0 or len(other) == 0:
            return np.array([], dtype=I64)
        n = max(len(self), len(other))
        for length in (len(self), len(other)):
            if length not in (1, n):
                raise ValueError("'x' and 'y' must have the same length, or one of them length 1")
        x = self if len(self) == n else self._take(np.zeros(n, dtype=np.int64))
        y = other if len(other) == n else other._take(np.zeros(n, dtype=np.int64))
        d = _nearest.distance_k(x._start, x._width, y._start, y._width)
        bad = x.seqnames != y.seqnames
        if not ignore_strand:
            bad = bad | (x._strand.astype(int) + y._strand.astype(int) == 1)
        return np.where(bad, NO_MATCH, d)

    # ----- coverage -----------------------------------------------------------------
    def coverage(
        self, shift: Any = 0, width: int | Mapping[str, int] | None = None, weight: Any = 1
    ) -> dict[str, Rle]:
        """R ``coverage``: ``{seqlevel: Rle}`` for every seqlevel, ignoring strand.
        Default length per level: its seqlength if known, else the largest end."""
        n = len(self)
        sh = V.recycle(V.as_int_array(shift, "shift"), n, "shift") if n else np.array([], dtype=I64)
        wt = np.asarray(weight)
        wt = np.repeat(wt, n) if wt.ndim == 0 else wt
        out: dict[str, Rle] = {}
        for i, name in enumerate(self._seqinfo.seqnames):
            rows = np.flatnonzero(self._seq == i)
            if isinstance(width, Mapping):
                length = width.get(name)
            elif width is not None:
                length = width
            else:
                length = self._seqinfo._lengths[i]
            out[name] = _coverage.coverage_k(
                self._start[rows], self._width[rows], sh[rows], length, wt[rows] if len(wt) else 1
            )
        return out

    # ----- ordering ----------------------------------------------------------
    def _order_keys(self, ignore_strand: bool = False) -> list[np.ndarray]:
        keys = [self._seq.astype(I64)]
        if not ignore_strand:
            keys.append(self._strand.astype(I64))
        return keys + [self._start, self._width]

    def order(self, decreasing: bool = False, ignore_strand: bool = False) -> np.ndarray:
        """R ``order``: by seqlevel (in seqinfo order), strand (``+``, ``-``,
        ``*``), start, then width; 0-based."""
        return _order.order_k(self._order_keys(ignore_strand), decreasing)

    def sort(self, decreasing: bool = False, ignore_strand: bool = False) -> GRanges:
        return self._take(self.order(decreasing, ignore_strand))

    def rank(self, ties_method: str = "first", ignore_strand: bool = False) -> np.ndarray:
        """R ``rank`` (1-based). Default ``"first"`` (R's generic has no default)."""
        return _order.rank_k(self._order_keys(ignore_strand), ties_method)

    def duplicated(self, from_last: bool = False, ignore_strand: bool = False) -> np.ndarray:
        """True for ranges equal in (seqname, strand, start, width) to an earlier
        one; names and mcols are ignored."""
        return _order.duplicated_k(self._order_keys(ignore_strand), from_last)

    def unique(self, from_last: bool = False, ignore_strand: bool = False) -> GRanges:
        return self._take(np.flatnonzero(~self.duplicated(from_last, ignore_strand)))

    # ----- comparison / display ----------------------------------------------
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GRanges):
            return NotImplemented
        if self._seqinfo != other._seqinfo:
            return False
        if not (
            np.array_equal(self._seq, other._seq)
            and np.array_equal(self._strand, other._strand)
            and np.array_equal(self._start, other._start)
            and np.array_equal(self._width, other._width)
        ):
            return False
        if (self._names is None) != (other._names is None):
            return False
        if self._names is not None and other._names is not None and list(self._names) != list(other._names):
            return False
        return self._mcols.equals(other._mcols) if (self._mcols.width or other._mcols.width) else True

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        n = len(self)
        rows = []
        sn, st = self.seqnames, self.strand
        for i in range(min(n, 10)):
            nm = f"{self._names[i]}: " if self._names is not None else ""
            rows.append(f"{nm}{sn[i]}:{self._start[i]}-{self._start[i] + self._width[i] - 1}:{st[i]}")
        more = ", ..." if n > 10 else ""
        extra = f", mcols={self._mcols.columns}" if self._mcols.width else ""
        return f"GRanges(n={n}{extra}, [{'; '.join(rows)}{more}], seqlevels={list(self.seqlevels)})"

    def to_tuples(self) -> list[tuple[str, int, int, str]]:
        """``[(seqname, start, end, strand), ...]`` as plain Python values (handy for tests)."""
        return list(zip(self.seqnames.tolist(), self._start.tolist(), self.end.tolist(), self.strand.tolist()))


_EMPTY_PAIRS = (np.array([], dtype=np.int64), np.array([], dtype=np.int64))


def _as_strand_codes(strand: str | Sequence[str] | None, n: int) -> np.ndarray:
    if strand is None:
        return np.full(n, STAR, dtype=np.int8)
    vals = [strand] if isinstance(strand, str) else [str(s) for s in strand]
    for v in vals:
        if v not in _STRAND_CODE:
            raise ValueError("strand values must be in '+' '-' '*'")
    if len(vals) == 1 and n != 1:
        vals = vals * n
    if len(vals) != n:
        raise ValueError(f"'strand' must have the length of the object ({n}) or length 1")
    return np.array([_STRAND_CODE[v] for v in vals], dtype=np.int8)


def _recode(codes: np.ndarray, old: Seqinfo, new: Seqinfo) -> np.ndarray:
    lookup = np.array([new.index(s) for s in old.seqnames], dtype=np.int32)
    return lookup[codes] if len(codes) else codes.astype(np.int32)


def _ro(a: np.ndarray) -> np.ndarray:
    v = a.view()
    v.flags.writeable = False
    return v
