"""``GRanges``: IRanges + seqnames + strand + seqinfo (spec 95-granges.md)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Self

import numpy as np
import numpy.typing as npt
import polars as pl

from ranges_e3 import _kernels, _mcols, _subscript
from ranges_e3._types import BoolArray, IntArray
from ranges_e3.hits import Hits
from ranges_e3.iranges import IRanges, per_range
from ranges_e3.rle import Rle
from ranges_e3.seqinfo import Seqinfo

STRANDS = ("+", "-", "*")  # code order = sort order (R factor order)
PLUS, MINUS, STAR = 0, 1, 2


def _strand_codes(strand: Any, n: int) -> npt.NDArray[np.int8]:
    """GR.1-2: scalar or per-range strand strings -> codes."""
    if strand is None:
        return np.full(n, STAR, dtype=np.int8)
    values = [strand] if isinstance(strand, str) else list(strand)
    if any(v not in STRANDS for v in values):
        raise ValueError("strand values must be '+', '-' or '*'")
    if len(values) not in (1, n):
        raise ValueError(f"strand has length {len(values)}, expected 1 or {n}")
    codes = np.array([STRANDS.index(v) for v in values], dtype=np.int8)
    return np.full(n, codes[0], dtype=np.int8) if len(values) == 1 else codes


class GRanges:
    """Immutable genomic ranges; seqnames are codes into ``seqinfo``'s levels."""

    __slots__ = ("_seq", "_ranges", "_strand", "_seqinfo", "_names", "_mcols")

    def __init__(
        self,
        seqnames: Sequence[str] | str,
        ranges: IRanges,
        strand: Sequence[str] | str | None = None,
        seqinfo: Seqinfo | None = None,
        names: Sequence[str | None] | None = None,
        mcols: pl.DataFrame | None = None,
    ) -> None:
        seq = [seqnames] if isinstance(seqnames, str) else list(seqnames)
        n = max(len(seq), len(ranges))
        if len(seq) not in (1, n) or len(ranges) not in (1, n):
            raise ValueError("seqnames and ranges lengths do not match")
        if seqinfo is None:
            seqinfo = Seqinfo(list(dict.fromkeys(seq)))
        unknown = sorted({s for s in seq if s not in seqinfo})
        if unknown:
            raise ValueError(f"seqnames not in seqinfo: {unknown[:3]}")
        codes = np.array([seqinfo.index(s) for s in seq], dtype=np.int32)
        if len(seq) == 1:
            codes = np.full(n, codes[0], dtype=np.int32)
        if len(ranges) == 1 and n != 1:
            ranges = ranges[np.zeros(n, dtype=np.int64)]
        if names is None:
            names = ranges.names
        self._seq: npt.NDArray[np.int32] = codes
        self._ranges = ranges.bare()
        self._strand = _strand_codes(strand, n)
        self._seqinfo = seqinfo
        self._names = IRanges._check_names(names, n)
        self._mcols = _mcols.validate(mcols, n)

    @classmethod
    def _new(
        cls,
        seq: npt.NDArray[np.int32],
        ranges: IRanges,
        strand: npt.NDArray[np.int8],
        seqinfo: Seqinfo,
        names: list[str | None] | None = None,
        mcols: pl.DataFrame | None = None,
    ) -> Self:
        """Trusted constructor (``ranges`` bare, codes valid)."""
        obj = cls.__new__(cls)
        obj._seq, obj._ranges, obj._strand = seq, ranges, strand
        obj._seqinfo, obj._names, obj._mcols = seqinfo, names, mcols
        return obj

    # --- accessors (GR.3, SEQI.3) ---

    def __len__(self) -> int:
        return len(self._seq)

    @property
    def seqnames(self) -> npt.NDArray[np.str_]:
        levels = np.array(self._seqinfo.seqnames + [""], dtype=np.str_)
        return levels[self._seq]

    @property
    def strand(self) -> npt.NDArray[np.str_]:
        return np.array(STRANDS, dtype=np.str_)[self._strand]

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
    def names(self) -> list[str | None] | None:
        return None if self._names is None else list(self._names)

    @property
    def mcols(self) -> pl.DataFrame:
        return _mcols.public(self._mcols)

    @property
    def seqinfo(self) -> Seqinfo:
        return self._seqinfo

    @property
    def seqlevels(self) -> list[str]:
        return self._seqinfo.seqnames

    @property
    def seqlengths(self) -> list[int | None]:
        return self._seqinfo.seqlengths

    @property
    def is_circular(self) -> list[bool | None]:
        return self._seqinfo.is_circular

    @property
    def genome(self) -> list[str | None]:
        return self._seqinfo.genome

    def with_names(self, names: Sequence[str | None] | None) -> Self:
        return self._rebuild(names=IRanges._check_names(names, len(self)))

    def with_mcols(self, mcols: pl.DataFrame | None) -> Self:
        return self._rebuild(mcols=_mcols.validate(mcols, len(self)))

    def _rebuild(self, **changes: Any) -> Self:
        fields = {
            "seq": self._seq, "ranges": self._ranges, "strand": self._strand,
            "seqinfo": self._seqinfo, "names": self._names, "mcols": self._mcols,
        }
        fields.update(changes)
        return self._new(**fields)

    # --- subsetting (GR.4) ---

    def __getitem__(self, index: Any) -> Self:
        return self._take(_subscript.resolve(index, len(self), self._names))

    def _take(self, pos: IntArray) -> Self:
        names = None if self._names is None else [self._names[i] for i in pos]
        return self._new(
            self._seq[pos], self._ranges[pos], self._strand[pos], self._seqinfo, names,
            _mcols.take(self._mcols, pos),
        )

    # --- intra-range (GR.7-8) ---

    def _like(self, start: IntArray, width: IntArray) -> Self:
        return self._rebuild(ranges=IRanges._new(start, width))

    def _mirrored(self, ignore_strand: bool) -> BoolArray:
        """Ranges whose start/end roles swap: ``-`` strand unless strands are ignored."""
        n = len(self)
        return np.zeros(n, dtype=np.bool_) if ignore_strand else self._strand == MINUS

    def shift(self, shift: Any) -> Self:
        """GR.7."""
        r = self._ranges
        return self._like(*_kernels.shift(r._start, r._width, per_range(shift, len(self), "shift")))

    def narrow(self, start: Any = None, end: Any = None, width: Any = None) -> Self:
        """GR.7."""
        r = self._ranges.narrow(start, end, width)
        return self._rebuild(ranges=r)

    def resize(self, width: Any, fix: str = "start", ignore_strand: bool = False) -> Self:
        """GR.8."""
        if fix not in ("start", "end", "center"):
            raise ValueError("fix must be 'start', 'end' or 'center'")
        n, r = len(self), self._ranges
        if fix == "center":
            codes = np.full(n, _kernels.FIX_CENTER, dtype=np.int64)
        else:
            plain, flipped = (
                (_kernels.FIX_START, _kernels.FIX_END) if fix == "start" else (_kernels.FIX_END, _kernels.FIX_START)
            )
            codes = np.where(self._mirrored(ignore_strand), flipped, plain).astype(np.int64)
        return self._like(*_kernels.resize(r._start, r._width, per_range(width, n, "width"), codes))

    def flank(self, width: Any, start: Any = True, both: bool = False, ignore_strand: bool = False) -> Self:
        """GR.8."""
        n, r = len(self), self._ranges
        at_start = np.broadcast_to(np.asarray(start, dtype=np.bool_), n) ^ self._mirrored(ignore_strand)
        return self._like(*_kernels.flank(r._start, r._width, per_range(width, n, "width"), at_start, both))

    def promoters(self, upstream: Any = 2000, downstream: Any = 200) -> Self:
        """GR.8."""
        n, r = len(self), self._ranges
        return self._like(
            *_kernels.promoters(
                r._start, r._width, per_range(upstream, n, "upstream"), per_range(downstream, n, "downstream"),
                self._mirrored(False),
            )
        )

    # --- grouping helpers ---

    def _group_keys(self, ignore_strand: bool) -> IntArray:
        """Sortable group id per range: ``(seqlevel, strand)`` or ``seqlevel`` alone (GR.10)."""
        strand = np.full(len(self), STAR, dtype=np.int64) if ignore_strand else self._strand.astype(np.int64)
        return self._seq.astype(np.int64) * 3 + strand

    @staticmethod
    def _groups(keys: IntArray) -> list[tuple[int, IntArray]]:
        """``[(key, indices)]`` ascending by key; indices ascending."""
        if len(keys) == 0:
            return []
        order = np.argsort(keys, kind="stable")
        sorted_keys = keys[order]
        cuts = np.flatnonzero(np.diff(sorted_keys)) + 1
        return [(int(keys[c[0]]), c) for c in np.split(order, cuts)]

    def _from_groups(
        self, parts: list[tuple[int, IntArray, IntArray]], info: Seqinfo | None = None
    ) -> Self:
        """Assemble ``(key, start, end)`` parts (already in output order) into a GRanges."""
        if parts:
            keys = np.concatenate([np.full(len(s), k, dtype=np.int64) for k, s, _ in parts])
            start = np.concatenate([s for _, s, _ in parts])
            end = np.concatenate([e for _, _, e in parts])
        else:
            keys, start, end = _mcols.empty_int(), _mcols.empty_int(), _mcols.empty_int()
        return self._new(
            (keys // 3).astype(np.int32), IRanges._new(start, end - start + 1), (keys % 3).astype(np.int8),
            info if info is not None else self._seqinfo,
        )

    # --- inter-range (GR.10-11) ---

    def range(self, ignore_strand: bool = False) -> Self:
        """GR.10."""
        r = self._ranges
        parts = []
        for key, idx in self._groups(self._group_keys(ignore_strand)):
            lo, hi = _kernels.span(r._start[idx], r._end[idx])  # type: ignore[misc]
            parts.append((key, np.array([lo], dtype=np.int64), np.array([hi], dtype=np.int64)))
        return self._from_groups(parts)

    def reduce(
        self, min_gapwidth: int = 1, drop_empty_ranges: bool = False, with_revmap: bool = False,
        ignore_strand: bool = False,
    ) -> Self:
        """GR.10."""
        r = self._ranges
        parts = []
        revmap: list[IntArray] = []
        for key, idx in self._groups(self._group_keys(ignore_strand)):
            gs_, ge, order, bounds = _kernels.reduce_scan(
                r._start[idx], r._end[idx], int(min_gapwidth), drop_empty_ranges
            )
            parts.append((key, gs_, ge))
            if with_revmap:
                groups = [idx[g] for g in _kernels.split_groups(order, bounds)]
                if drop_empty_ranges:
                    groups = [g[r._width[g] > 0] for g in groups]
                revmap.extend(groups)
        out = self._from_groups(parts)
        return out.with_mcols(pl.DataFrame([_mcols.revmap_series(revmap)])) if with_revmap else out

    def disjoin(self, with_revmap: bool = False, ignore_strand: bool = False) -> Self:
        """GR.10."""
        r = self._ranges
        parts = []
        revmap: list[IntArray] = []
        for key, idx in self._groups(self._group_keys(ignore_strand)):
            ps, pe, rm = _kernels.disjoin(r._start[idx], r._end[idx], with_revmap)
            parts.append((key, ps, pe))
            revmap.extend(idx[g] for g in rm)
        out = self._from_groups(parts)
        return out.with_mcols(pl.DataFrame([_mcols.revmap_series(revmap)])) if with_revmap else out

    def gaps(self, start: int | None = None, end: int | None = None, ignore_strand: bool = False) -> Self:
        """GR.11."""
        r = self._ranges
        strands = (STAR,) if ignore_strand else (PLUS, MINUS, STAR)
        keys = self._group_keys(ignore_strand)
        parts = []
        for level, length in enumerate(self.seqlengths):
            for strand in strands:
                key = level * 3 + strand
                idx = np.flatnonzero(keys == key)
                lo = 1 if start is None else int(start)
                if end is not None:
                    hi: int | None = int(end)
                elif length is not None:
                    hi = length
                else:
                    positive = idx[r._width[idx] > 0]
                    hi = int(r._end[positive].max()) if len(positive) else None
                if hi is None:
                    continue
                gs_, ge = _kernels.gaps(r._start[idx], r._end[idx], lo, hi)
                if len(gs_):
                    parts.append((key, gs_, ge))
        return self._from_groups(parts)

    # --- set operations (GR.12-13) ---

    def _aligned(self, y: GRanges) -> tuple[Seqinfo, npt.NDArray[np.int32]]:
        """Merged seqinfo and ``y``'s seqname codes expressed in it."""
        info = self._seqinfo.merge(y._seqinfo)
        remap = np.array([info.index(name) for name in y._seqinfo.seqnames] + [0], dtype=np.int32)
        return info, remap[y._seq]

    def _set_op(self, y: GRanges, op: str, ignore_strand: bool) -> Self:
        info, y_seq = self._aligned(y)
        strand_y = np.full(len(y), STAR, dtype=np.int64) if ignore_strand else y._strand.astype(np.int64)
        kx = self._group_keys(ignore_strand)
        ky = y_seq.astype(np.int64) * 3 + strand_y
        gx, gy = dict(self._groups(kx)), dict(self._groups(ky))
        wanted = {"union": sorted(set(gx) | set(gy)), "intersect": sorted(set(gx) & set(gy)), "setdiff": sorted(gx)}[op]
        e = _mcols.empty_int()
        parts = []
        for key in wanted:
            xi, yi = gx.get(key, e), gy.get(key, e)
            s, en = _kernels.set_operation(
                self._ranges._start[xi], self._ranges._end[xi], y._ranges._start[yi], y._ranges._end[yi], op
            )
            if op != "union":  # R clips intersect/setdiff results to [1, seqlength] (GR.12)
                limit = info.seqlengths[key // 3]
                if limit is not None:
                    s, en = np.maximum(s, 1), np.minimum(en, limit)
                    ok = s <= en
                    s, en = s[ok], en[ok]
            if len(s):
                parts.append((key, s, en))
        return self._from_groups(parts, info)

    def union(self, y: GRanges, ignore_strand: bool = False) -> Self:
        """GR.12."""
        return self._set_op(y, "union", ignore_strand)

    def intersect(self, y: GRanges, ignore_strand: bool = False) -> Self:
        """GR.12."""
        return self._set_op(y, "intersect", ignore_strand)

    def setdiff(self, y: GRanges, ignore_strand: bool = False) -> Self:
        """GR.12."""
        return self._set_op(y, "setdiff", ignore_strand)

    def _comparable(self, y: GRanges, y_seq: npt.NDArray[np.int32], ignore_strand: bool) -> BoolArray:
        """GR.6 for aligned pairs (``y_seq`` in this object's seqinfo codes)."""
        same = self._seq == y_seq
        if ignore_strand:
            return same
        return same & ((self._strand == y._strand) | (self._strand == STAR) | (y._strand == STAR))

    def _parallel(
        self, y: GRanges, ignore_strand: bool, kernel: Any, fallback: str, y_strand: Any
    ) -> Self:
        """Pairwise op on comparable pairs; ``y_strand(ok, xs, xe, ys, ye)`` flags pairs that take ``y``'s strand."""
        if len(self) != len(y):
            raise ValueError("parallel operations need equal lengths (SET.2)")
        info, y_seq = self._aligned(y)
        ok = self._comparable(y, y_seq, ignore_strand)
        xs, xe = self._ranges._start, self._ranges._end
        ys, ye = y._ranges._start, y._ranges._end
        s, e = xs.copy(), xe.copy()
        if ok.any():
            s[ok], e[ok] = kernel(xs[ok], xe[ok], ys[ok], ye[ok])
        if not ok.all():
            if fallback == "refuse":
                raise ValueError("pair of ranges is not comparable (different seqname or strand)")
            if fallback == "empty":
                e[~ok] = xs[~ok] - 1
        strand = self._strand.copy()
        if not ignore_strand:
            take_y = (self._strand == STAR) & y_strand(ok, xs, xe, ys, ye)
            strand = np.where(take_y, y._strand, self._strand).astype(np.int8)
        return self._new(self._seq, IRanges._new(s, e - s + 1), strand, info, self._names)

    def punion(self, y: GRanges, fill_gap: bool = False, ignore_strand: bool = False) -> Self:
        """GR.13."""
        return self._parallel(
            y, ignore_strand, lambda a, b, c, d: _kernels.punion(a, b, c, d, fill_gap), "refuse",
            lambda ok, *_: np.zeros(len(ok), dtype=np.bool_),
        )

    def pintersect(self, y: GRanges, ignore_strand: bool = False) -> Self:
        """GR.13 (pairs without overlap give a zero-width range at ``x``'s start)."""
        return self._parallel(
            y, ignore_strand, lambda a, b, c, d: _kernels.pintersect(a, b, c, d, "start.x"), "empty",
            lambda ok, xs, xe, ys, ye: ok & (np.minimum(xe, ye) >= np.maximum(xs, ys) - 1),
        )

    def psetdiff(self, y: GRanges, ignore_strand: bool = False) -> Self:
        """GR.13."""
        return self._parallel(
            y, ignore_strand, _kernels.psetdiff, "keep", lambda ok, *_: np.ones(len(ok), dtype=np.bool_)
        )

    # --- coverage (GR.17) ---

    def coverage(self, shift: Any = 0, width: int | None = None) -> dict[str, Rle]:
        """GR.17: per-level run-length coverage, strand ignored."""
        r = self._ranges
        shifts = per_range(shift, len(self), "shift")
        out: dict[str, Rle] = {}
        for level, (name, length) in enumerate(zip(self.seqlevels, self.seqlengths)):
            idx = np.flatnonzero(self._seq == level)
            size = width if width is not None else length
            lengths, values = _kernels.coverage_runs(r._start[idx], r._end[idx], shifts[idx], size)
            out[name] = Rle(lengths, values)
        return out

    # --- overlaps (GR.14) ---

    def _blocks(
        self, y: GRanges, ignore_strand: bool, split_star: bool = False
    ) -> list[tuple[IntArray, IntArray, bool]]:
        """Comparable (query indices, subject indices, mirrored) blocks (GR.6).

        Each query is in one block, except with ``split_star`` where a ``*`` query is in two:
        against ``+``/``*`` subjects (not mirrored) and against ``-`` subjects (mirrored, GR.15).
        """
        _, y_seq = self._aligned(y)
        blocks: list[tuple[IntArray, IntArray, bool]] = []
        for seq in np.unique(self._seq):
            q_in = np.flatnonzero(self._seq == seq)
            s_in = np.flatnonzero(y_seq == seq)
            if len(s_in) == 0:
                continue
            if ignore_strand:
                blocks.append((q_in, s_in, False))
                continue
            s_strand = y._strand[s_in]
            for strand in (PLUS, MINUS, STAR):
                q_idx = q_in[self._strand[q_in] == strand]
                if len(q_idx) == 0:
                    continue
                if strand == STAR and split_star:
                    blocks.append((q_idx, s_in[s_strand != MINUS], False))
                    blocks.append((q_idx, s_in[s_strand == MINUS], True))
                    continue
                keep = np.ones(len(s_in), dtype=np.bool_) if strand == STAR else (
                    (s_strand == strand) | (s_strand == STAR)
                )
                blocks.append((q_idx, s_in[keep], strand == MINUS))
        return blocks

    def _block_pairs(
        self, y: GRanges, ignore_strand: bool, kernel: Any, split_star: bool = False
    ) -> tuple[IntArray, IntArray]:
        """Run ``kernel(q_start, q_end, s_start, s_end, mirrored, q_idx, s_idx)`` per block; sort pairs."""
        qs, qe = self._ranges._start, self._ranges._end
        ss, se = y._ranges._start, y._ranges._end
        out_q: list[IntArray] = [_mcols.empty_int()]
        out_s: list[IntArray] = [_mcols.empty_int()]
        for q_idx, s_idx, mirrored in self._blocks(y, ignore_strand, split_star):
            qi, si = kernel(qs[q_idx], qe[q_idx], ss[s_idx], se[s_idx], mirrored, q_idx, s_idx)
            out_q.append(q_idx[qi])
            out_s.append(s_idx[si])
        qi, si = np.concatenate(out_q), np.concatenate(out_s)
        order = np.lexsort((si, qi))
        return qi[order], si[order]

    def _overlap_pairs(
        self, subject: GRanges, max_gap: int, min_overlap: int, type: str, ignore_strand: bool
    ) -> tuple[IntArray, IntArray]:
        _kernels.check_overlap_args(int(max_gap), int(min_overlap), type)
        return self._block_pairs(
            subject, ignore_strand,
            lambda a, b, c, d, _m, _q, _s: _kernels.find_overlaps(a, b, c, d, int(max_gap), int(min_overlap), type),
        )

    def find_overlaps(
        self, subject: GRanges, max_gap: int = -1, min_overlap: int = 0, type: str = "any",
        select: str = "all", ignore_strand: bool = False,
    ) -> Hits | np.ma.MaskedArray:  # type: ignore[type-arg]
        """GR.14."""
        qi, si = self._overlap_pairs(subject, max_gap, min_overlap, type, ignore_strand)
        if select == "all":
            return Hits(qi, si, len(self), len(subject))
        return _kernels.select_pairs(qi, si, len(self), select)

    def count_overlaps(
        self, subject: GRanges, max_gap: int = -1, min_overlap: int = 0, type: str = "any",
        ignore_strand: bool = False,
    ) -> IntArray:
        """GR.14."""
        qi, _ = self._overlap_pairs(subject, max_gap, min_overlap, type, ignore_strand)
        return np.bincount(qi, minlength=len(self)).astype(np.int64)

    def overlaps_any(
        self, subject: GRanges, max_gap: int = -1, min_overlap: int = 0, type: str = "any",
        ignore_strand: bool = False,
    ) -> BoolArray:
        """GR.14."""
        return self.count_overlaps(subject, max_gap, min_overlap, type, ignore_strand) > 0

    def subset_by_overlaps(
        self, subject: GRanges, max_gap: int = -1, min_overlap: int = 0, type: str = "any",
        invert: bool = False, ignore_strand: bool = False,
    ) -> Self:
        """GR.14."""
        hit = self.overlaps_any(subject, max_gap, min_overlap, type, ignore_strand)
        return self._take(np.flatnonzero(~hit if invert else hit).astype(np.int64))

    # --- nearest / precede / follow / distance (GR.15) ---

    def _neighbors(
        self, subject: GRanges | None, select: str, allowed: tuple[str, ...], ignore_strand: bool, kernel: Any,
        pick: str, directional: bool = False,
    ) -> Hits | np.ma.MaskedArray:  # type: ignore[type-arg]
        if select not in allowed:
            raise ValueError(f"select must be one of {allowed}")
        target = self if subject is None else subject
        qi, si = self._block_pairs(target, ignore_strand, kernel, split_star=directional)
        if directional:
            qi, si = self._closest_per_query(target, qi, si)
        if select == "all":
            return Hits(qi, si, len(self), len(target))
        return _kernels.select_pairs(qi, si, len(self), pick)

    def _closest_per_query(self, target: GRanges, qi: IntArray, si: IntArray) -> tuple[IntArray, IntArray]:
        """Keep, per query, the pairs at minimal distance (merges the two ``*``-query blocks)."""
        if len(qi) == 0:
            return qi, si
        r, t = self._ranges, target._ranges
        d = _kernels.distance(r._start[qi], r._end[qi], t._start[si], t._end[si])
        starts = np.flatnonzero(np.concatenate([[True], qi[1:] != qi[:-1]]))
        counts = np.diff(np.append(starts, len(qi)))
        best = np.repeat(np.minimum.reduceat(d, starts), counts)
        keep = d == best
        return qi[keep], si[keep]

    def nearest(
        self, subject: GRanges | None = None, select: str = "arbitrary", ignore_strand: bool = False
    ) -> Hits | np.ma.MaskedArray:  # type: ignore[type-arg]
        """GR.15."""
        def kernel(a: IntArray, b: IntArray, c: IntArray, d: IntArray, _m: bool, q_idx: IntArray, s_idx: IntArray) -> Any:
            ids = (q_idx, s_idx) if subject is None else None
            return _kernels.nearest_pairs(a, b, c, d, ids)

        return self._neighbors(subject, select, ("arbitrary", "all"), ignore_strand, kernel, "first")

    def precede(
        self, subject: GRanges | None = None, select: str = "first", ignore_strand: bool = False
    ) -> Hits | np.ma.MaskedArray:  # type: ignore[type-arg]
        """GR.15 (a ``-`` query looks toward lower coordinates)."""
        def kernel(a: IntArray, b: IntArray, c: IntArray, d: IntArray, mirrored: bool, *_: Any) -> Any:
            return _kernels.precede_pairs(-a, -d) if mirrored else _kernels.precede_pairs(b, c)

        return self._neighbors(subject, select, ("first", "all"), ignore_strand, kernel, "first", directional=True)

    def follow(
        self, subject: GRanges | None = None, select: str = "first", ignore_strand: bool = False
    ) -> Hits | np.ma.MaskedArray:  # type: ignore[type-arg]
        """GR.15 (``first`` = highest index of a tie, ADR-0010)."""
        def kernel(a: IntArray, b: IntArray, c: IntArray, d: IntArray, mirrored: bool, *_: Any) -> Any:
            return _kernels.follow_pairs(-b, -c) if mirrored else _kernels.follow_pairs(a, d)

        if select not in ("first", "all"):
            raise ValueError("select must be 'first' or 'all'")
        return self._neighbors(
            subject, select, ("first", "all"), ignore_strand, kernel, "last", directional=True
        )

    def distance(self, y: GRanges, ignore_strand: bool = False) -> np.ma.MaskedArray:  # type: ignore[type-arg]
        """GR.15: masked (NA) where the pair is not comparable."""
        n, m = len(self), len(y)
        if n != m and 1 not in (n, m):
            raise ValueError("distance needs equal lengths or one of length 1")
        if 0 in (n, m):
            return _kernels.masked_ints(_mcols.empty_int(), np.zeros(0, dtype=np.bool_))
        size = max(n, m)
        xi = np.zeros(size, dtype=np.int64) if n == 1 else np.arange(size)
        yi = np.zeros(size, dtype=np.int64) if m == 1 else np.arange(size)
        x, yy = self[xi], y[yi]
        _, y_seq = x._aligned(yy)
        ok = x._comparable(yy, y_seq, ignore_strand)
        d = _kernels.distance(x._ranges._start, x._ranges._end, yy._ranges._start, yy._ranges._end)
        return _kernels.masked_ints(d, ~ok)

    # --- ordering (GR.16) ---

    def _keys(self) -> list[IntArray]:
        return [self._seq.astype(np.int64), self._strand.astype(np.int64), self._ranges._start, self._ranges._end]

    def order(self, decreasing: bool = False) -> IntArray:
        """GR.16."""
        return _kernels.order_by(self._keys(), decreasing)

    def sort(self, decreasing: bool = False) -> Self:
        """GR.16."""
        return self._take(self.order(decreasing))

    def rank(self, ties_method: str = "average") -> IntArray | npt.NDArray[np.float64]:
        """GR.16."""
        return _kernels.rank_by(self._keys(), ties_method)

    def duplicated(self, from_last: bool = False) -> BoolArray:
        """GR.16."""
        return _kernels.duplicated_by(self._keys(), from_last)

    def unique(self) -> Self:
        """GR.16."""
        return self._take(np.flatnonzero(~self.duplicated()).astype(np.int64))

    # --- trim (GR.9) ---

    def trim(self) -> Self:
        """Clip ranges to ``[1, seqlength]`` on non-circular levels of known length."""
        lengths = np.array(
            [-1 if (L is None or c) else L for L, c in zip(self.seqlengths, self.is_circular)],
            dtype=np.int64,
        )
        limit = lengths[self._seq]
        active = limit >= 0
        s, e = self._ranges._start, self._ranges.end
        new_s = np.where(active, np.minimum(np.maximum(s, 1), limit + 1), s)
        new_e = np.where(active, np.maximum(np.minimum(e, limit), 0), e)
        return self._rebuild(ranges=IRanges._new(new_s, new_e - new_s + 1))

    def __repr__(self) -> str:
        shown = ", ".join(
            f"{sn}:{s}-{e}:{st}"
            for sn, s, e, st in zip(self.seqnames[:4], self.start[:4], self.end[:4], self.strand[:4])
        )
        return f"GRanges(n={len(self)}: {shown}{', ...' if len(self) > 4 else ''})"


def concat_granges(objs: Sequence[GRanges]) -> GRanges:
    """GR.5: append, merging seqinfo left to right and remapping seqname codes."""
    if not objs:
        raise ValueError("concat needs at least one object")
    info = objs[0].seqinfo
    for o in objs[1:]:
        info = info.merge(o.seqinfo)
    codes = []
    for o in objs:
        remap = np.array([info.index(name) for name in o.seqinfo.seqnames] + [0], dtype=np.int32)
        codes.append(remap[o._seq])
    names: list[str | None] | None = None
    if any(o._names is not None for o in objs):
        names = []
        for o in objs:
            names.extend(o._names if o._names is not None else [""] * len(o))
    return GRanges._new(
        np.concatenate(codes).astype(np.int32),
        IRanges._new(
            np.concatenate([o._ranges._start for o in objs]),
            np.concatenate([o._ranges._width for o in objs]),
        ),
        np.concatenate([o._strand for o in objs]),
        info,
        names,
        _mcols.concat([(o._mcols, len(o)) for o in objs]),
    )
