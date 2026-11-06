"""GRanges: IRanges + seqnames + strand + seqinfo. Strand-aware where R is strand-aware."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import polars as pl

from . import _coverage, _inter, _mcols, _nearest, _ordering, _overlaps
from ._util import IntArray, as_bool_array, as_int64, positions_from_subscript, rep_len
from .hits import Hits
from .iranges import IRanges, _revmap_frame, _select_hits, _ties_result
from .rle import Factor, Rle
from .seqinfo import Seqinfo

STRAND_LEVELS = ("+", "-", "*")
PLUS, MINUS, STAR = 0, 1, 2
NA_INDEX = -1  # sentinel for "no hit" / NA in index and distance outputs

SeqnamesArg = Sequence[str] | Factor | None


def _strand_codes(strand: Any, n: int) -> npt.NDArray[np.int8]:
    if strand is None:
        return np.full(n, STAR, dtype=np.int8)
    if isinstance(strand, str):
        strand = [strand]
    vals = list(strand)
    if len(vals) == 1 and n != 1:
        vals = vals * n
    if len(vals) != n:
        raise ValueError(f"strand must have length {n} or length 1, got {len(vals)}")
    try:
        return np.array([STAR if v is None else STRAND_LEVELS.index(v) for v in vals], dtype=np.int8)
    except ValueError:
        raise ValueError("strand values must be in '+' '-' '*'") from None


def _compat(qstrand: npt.NDArray[np.int8] | int, sstrand: npt.NDArray[np.int8]) -> npt.NDArray[np.bool_]:
    """Strands are compatible when equal or either is '*'."""
    return (np.asarray(qstrand) == STAR) | (sstrand == STAR) | (np.asarray(qstrand) == sstrand)


class GRanges:
    """Immutable genomic ranges. `seqnames` are codes into `seqinfo`, `strand` 0 '+', 1 '-', 2 '*'."""

    __slots__ = ("_seq", "_strand", "_ranges", "_seqinfo")

    def __init__(
        self,
        seqnames: SeqnamesArg = None,
        ranges: IRanges | None = None,
        strand: Sequence[str | None] | str | None = None,
        *,
        seqinfo: Seqinfo | None = None,
        seqlengths: Mapping[str, int | None] | None = None,
        mcols: pl.DataFrame | None = None,
        **columns: Any,
    ) -> None:
        rng = ranges if ranges is not None else IRanges()
        n = len(rng)
        names, default_levels = self._resolve_seqnames(seqnames, n)
        if seqinfo is None:
            if seqlengths:
                seqinfo = Seqinfo(list(seqlengths), list(seqlengths.values()))
            else:
                seqinfo = Seqinfo(default_levels)
        index = {lv: i for i, lv in enumerate(seqinfo.seqnames)}
        if any(nm not in index for nm in names):
            raise ValueError("seqnames contains sequence names with no entries in seqinfo")
        codes = np.array([index[nm] for nm in names], dtype=np.int32)
        extra = mcols
        if columns:
            frame = pl.DataFrame({k: _series_values(v) for k, v in columns.items()})
            extra = frame if extra is None else pl.concat([extra, frame], how="horizontal")
        if extra is not None:
            rng = rng.with_mcols(extra)
        object.__setattr__(self, "_seq", codes)
        object.__setattr__(self, "_strand", _strand_codes(strand, n))
        object.__setattr__(self, "_ranges", rng)
        object.__setattr__(self, "_seqinfo", seqinfo)

    @staticmethod
    def _resolve_seqnames(seqnames: SeqnamesArg, n: int) -> tuple[list[str], list[str]]:
        """Return (per-range names, default levels): factor levels, else first-appearance order."""
        if seqnames is None:
            if n:
                raise ValueError("seqnames must be supplied")
            return [], []
        if isinstance(seqnames, Factor):
            names = seqnames.to_list()
            levels = list(seqnames.levels)
        else:
            names = list(seqnames)
            levels = list(dict.fromkeys(names))
        if any(v is None for v in names):
            raise ValueError("seqnames cannot contain NAs")
        if len(names) == 1 and n != 1:
            names = names * n
        if len(names) != n:
            raise ValueError(f"seqnames must have the length of the object to construct ({n}) or length 1")
        return names, levels  # type: ignore[return-value]

    def __setattr__(self, key: str, value: Any) -> None:
        raise AttributeError("GRanges is immutable")

    @classmethod
    def _make(
        cls,
        seq: npt.NDArray[np.integer[Any]],
        strand: npt.NDArray[np.integer[Any]],
        ranges: IRanges,
        seqinfo: Seqinfo,
    ) -> "GRanges":
        obj = object.__new__(cls)
        object.__setattr__(obj, "_seq", np.asarray(seq, dtype=np.int32))
        object.__setattr__(obj, "_strand", np.asarray(strand, dtype=np.int8))
        object.__setattr__(obj, "_ranges", ranges)
        object.__setattr__(obj, "_seqinfo", seqinfo)
        return obj

    # ---- accessors ---------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._ranges)

    @property
    def seqnames(self) -> Rle:
        return Rle.from_factor(Factor(self._seqinfo.seqnames, self._seq.astype(np.int64)))

    @property
    def strand(self) -> Rle:
        return Rle.from_factor(Factor(STRAND_LEVELS, self._strand.astype(np.int64)))

    @property
    def ranges(self) -> IRanges:
        """The ranges with names but without metadata columns."""
        return IRanges._make(self._ranges._s(), self._ranges._w(), self._ranges._names)

    @property
    def start(self) -> npt.NDArray[np.int32]:
        return self._ranges.start

    @property
    def end(self) -> npt.NDArray[np.int32]:
        return self._ranges.end

    @property
    def width(self) -> npt.NDArray[np.int32]:
        return self._ranges.width

    @property
    def names(self) -> list[str | None] | None:
        return self._ranges.names

    @property
    def mcols(self) -> pl.DataFrame:
        return self._ranges.mcols

    @property
    def has_mcols(self) -> bool:
        return self._ranges.has_mcols

    @property
    def seqinfo(self) -> Seqinfo:
        return self._seqinfo

    @property
    def seqlevels(self) -> list[str]:
        return list(self._seqinfo.seqnames)

    def seqname_strings(self) -> list[str]:
        return [self._seqinfo.seqnames[c] for c in self._seq.tolist()]

    def strand_strings(self) -> list[str]:
        return [STRAND_LEVELS[c] for c in self._strand.tolist()]

    def _s(self) -> IntArray:
        return self._ranges._s()

    def _w(self) -> IntArray:
        return self._ranges._w()

    def _with(
        self,
        *,
        ranges: IRanges | None = None,
        seq: Any = None,
        strand: Any = None,
        seqinfo: Seqinfo | None = None,
    ) -> "GRanges":
        return GRanges._make(
            self._seq if seq is None else seq,
            self._strand if strand is None else strand,
            self._ranges if ranges is None else ranges,
            self._seqinfo if seqinfo is None else seqinfo,
        )

    def with_names(self, names: Sequence[str | None] | None) -> "GRanges":
        return self._with(ranges=self._ranges.with_names(names))

    def with_mcols(self, mcols: pl.DataFrame | None) -> "GRanges":
        return self._with(ranges=self._ranges.with_mcols(mcols))

    # ---- seqinfo edits -------------------------------------------------------------

    def with_seqinfo(self, seqinfo: Seqinfo) -> "GRanges":
        """R's `seqinfo<-` with identical level sets: replaces the metadata, relabels data by name."""
        if set(seqinfo.seqnames) != set(self._seqinfo.seqnames):
            raise ValueError("the seqlevels in the supplied Seqinfo are not the same as in x")
        remap = np.array([seqinfo.seqnames.index(nm) for nm in self._seqinfo.seqnames], dtype=np.int32)
        return self._with(seq=remap[self._seq] if len(self) else self._seq, seqinfo=seqinfo)

    def with_seqlevels(self, value: Sequence[str]) -> "GRanges":
        """R's `seqlevels<-`: rename (nothing in common, same length), reorder, add, or drop unused levels."""
        new, mapping = self._seqinfo.relevel(value)
        remap = np.array([-1 if m is None else m for m in mapping], dtype=np.int32)
        if len(self):
            new_codes = remap[self._seq]
            if np.any(new_codes < 0):
                used = sorted({self._seqinfo.seqnames[c] for c in self._seq[new_codes < 0].tolist()})
                raise ValueError(f"seqlevels to be dropped are in use: {', '.join(used)}")
        else:
            new_codes = self._seq
        return self._with(seq=new_codes, seqinfo=new)

    def with_seqlengths(self, value: Any) -> "GRanges":
        return self._with(seqinfo=self._seqinfo.with_field("seqlengths", value))

    def with_is_circular(self, value: Any) -> "GRanges":
        return self._with(seqinfo=self._seqinfo.with_field("is_circular", value))

    def with_genome(self, value: Any) -> "GRanges":
        return self._with(seqinfo=self._seqinfo.with_field("genome", value))

    # ---- subsetting / combining -----------------------------------------------------

    def __getitem__(self, i: Any) -> "GRanges":
        idx = positions_from_subscript(i, len(self), self._ranges._names)
        return self._take(idx)

    def _take(self, idx: IntArray) -> "GRanges":
        return self._with(seq=self._seq[idx], strand=self._strand[idx], ranges=self._ranges._take(idx))

    def concat(self, *others: "GRanges") -> "GRanges":
        parts = [self, *others]
        if any(not isinstance(p, GRanges) for p in parts):
            raise TypeError("can only concatenate GRanges with GRanges")
        si = parts[0]._seqinfo
        for p in parts[1:]:
            si = si.merge(p._seqinfo)
        codes = np.concatenate([p._recode(si) for p in parts])
        strand = np.concatenate([p._strand for p in parts])
        ranges = parts[0]._ranges.concat(*[p._ranges for p in others])
        return GRanges._make(codes, strand, ranges, si)

    def _recode(self, si: Seqinfo) -> npt.NDArray[np.int32]:
        """Seqname codes of self expressed against the (superset) seqinfo `si`."""
        if si.seqnames[: len(self._seqinfo)] == self._seqinfo.seqnames:
            return self._seq
        remap = np.array([si.seqnames.index(nm) for nm in self._seqinfo.seqnames], dtype=np.int32)
        return remap[self._seq] if len(self) else self._seq

    # ---- grouping helper ---------------------------------------------------------------

    def _group_keys(self, ignore_strand: bool) -> IntArray:
        st = np.full(len(self), STAR, dtype=np.int64) if ignore_strand else self._strand.astype(np.int64)
        return self._seq.astype(np.int64) * 3 + st

    def _groups(self, ignore_strand: bool = False) -> list[tuple[int, int, IntArray]]:
        """(seq code, strand code, member indices) for groups present, ordered seqlevel then + - *."""
        keys = self._group_keys(ignore_strand)
        out: list[tuple[int, int, IntArray]] = []
        if len(keys) == 0:
            return out
        order = np.argsort(keys, kind="stable")
        sorted_keys = keys[order]
        cuts = np.nonzero(np.diff(sorted_keys))[0] + 1
        for chunk in np.split(order, cuts):
            k = int(keys[chunk[0]])
            out.append((k // 3, k % 3, chunk.astype(np.int64)))
        return out

    def _assemble(
        self,
        pieces: list[tuple[int, int, IntArray, IntArray]],
        seqinfo: Seqinfo | None = None,
        mcols: pl.DataFrame | None = None,
    ) -> "GRanges":
        """Build a GRanges from (seq, strand, start, width) group results, already in output order."""
        if pieces:
            seq = np.concatenate([np.full(len(p[2]), p[0]) for p in pieces])
            st = np.concatenate([np.full(len(p[2]), p[1]) for p in pieces])
            s = np.concatenate([p[2] for p in pieces])
            w = np.concatenate([p[3] for p in pieces])
        else:
            seq = st = s = w = np.zeros(0, dtype=np.int64)
        return GRanges._make(seq, st, IRanges._make(s, w, None, mcols), seqinfo or self._seqinfo)

    # ---- intra-range ----------------------------------------------------------------------

    def shift(self, shift: Any) -> "GRanges":
        return self._with(ranges=self._ranges.shift(shift))

    def narrow(self, start: Any = None, end: Any = None, width: Any = None) -> "GRanges":
        return self._with(ranges=self._ranges.narrow(start, end, width))

    def _minus(self, ignore_strand: bool) -> npt.NDArray[np.bool_]:
        return np.zeros(len(self), dtype=bool) if ignore_strand else self._strand == MINUS

    def resize(
        self, width: Any, fix: str | Sequence[str] = "start", use_names: bool = True, ignore_strand: bool = False
    ) -> "GRanges":
        fx = np.asarray([fix] if isinstance(fix, str) else list(fix))
        fx = rep_len(fx, len(self), "fix") if len(self) else fx
        minus = self._minus(ignore_strand)
        swapped = np.where(fx == "start", "end", np.where(fx == "end", "start", fx))
        eff = np.where(minus, swapped, fx) if len(self) else fx
        return self._with(ranges=self._ranges.resize(width, eff, use_names))

    def flank(
        self, width: Any, start: bool | Sequence[bool] = True, both: bool = False, ignore_strand: bool = False
    ) -> "GRanges":
        st = as_bool_array([start] if isinstance(start, bool) else start, len(self), "start")
        st = np.where(self._minus(ignore_strand), ~st, st)
        return self._with(ranges=self._ranges.flank(width, st, both))

    def promoters(self, upstream: Any = 2000, downstream: Any = 200) -> "GRanges":
        from . import _intra

        ns, nw = _intra.promoter_ranges(
            self._s(),
            self._w(),
            as_int64(upstream, "upstream"),
            as_int64(downstream, "downstream"),
            self._strand == MINUS,
        )
        return self._with(ranges=IRanges._make(ns, nw, self._ranges._names, self._ranges._mcols))

    def trim(self) -> "GRanges":
        """Clip ranges to 1..seqlength on non-circular sequences with known length."""
        if len(self) == 0:
            return self
        lens = np.array(
            [-1 if (v is None or c) else v for v, c in zip(self._seqinfo.seqlengths, self._seqinfo.is_circular)],
            dtype=np.int64,
        )[self._seq]
        known = lens >= 0
        s, e = self._s(), self._s() + self._w() - 1
        new_s = np.where(known, np.clip(s, 1, lens + 1), s)
        new_e = np.where(known, np.clip(e, 0, lens), e)
        new_e = np.maximum(new_e, new_s - 1)
        return self._with(ranges=IRanges._make(new_s, new_e - new_s + 1, self._ranges._names, self._ranges._mcols))

    # ---- inter-range ---------------------------------------------------------------------------

    def range(self, ignore_strand: bool = False) -> "GRanges":
        pieces = []
        for seq, st, idx in self._groups(ignore_strand):
            s, w = _inter.range_of(self._s()[idx], self._w()[idx])
            pieces.append((seq, st, s, w))
        return self._assemble(pieces)

    def reduce(
        self,
        drop_empty_ranges: bool = False,
        min_gapwidth: int = 1,
        with_revmap: bool = False,
        ignore_strand: bool = False,
    ) -> "GRanges":
        pieces = []
        revmaps: list[IntArray] = []
        for seq, st, idx in self._groups(ignore_strand):
            s, w, rev = _inter.reduce_ranges(self._s()[idx], self._w()[idx], min_gapwidth, drop_empty_ranges, with_revmap)
            pieces.append((seq, st, s, w))
            if rev is not None:
                revmaps.extend(idx[r] for r in rev)
        return self._assemble(pieces, mcols=_revmap_frame(revmaps) if with_revmap else None)

    def disjoin(self, with_revmap: bool = False, ignore_strand: bool = False) -> "GRanges":
        pieces = []
        revmaps: list[IntArray] = []
        for seq, st, idx in self._groups(ignore_strand):
            s, w, rev = _inter.disjoin_ranges(self._s()[idx], self._w()[idx], with_revmap)
            pieces.append((seq, st, s, w))
            if rev is not None:
                revmaps.extend(idx[r] for r in rev)
        return self._assemble(pieces, mcols=_revmap_frame(revmaps) if with_revmap else None)

    def gaps(
        self,
        start: int | Sequence[int] | Mapping[str, int] | None = None,
        end: int | Sequence[int] | Mapping[str, int | None] | None = None,
        ignore_strand: bool = False,
    ) -> "GRanges":
        """Gaps for every seqlevel and strand (or '*' only with ignore_strand), within [start, end].

        start defaults to 1, end to the seqlength; with unknown seqlength the window ends at the
        last range of the group (a group without ranges then yields nothing).
        """
        si = self._seqinfo
        nlev = len(si)

        def per_level(arg: Any, default: list[Any]) -> list[Any]:
            if arg is None:
                return default
            if isinstance(arg, Mapping):
                return [arg.get(nm, default[i]) for i, nm in enumerate(si.seqnames)]
            vals = np.asarray(arg, dtype=object).reshape(-1)
            return [None if v is None else int(v) for v in rep_len(vals, nlev, "start/end")] if nlev else []

        starts = per_level(start, [1] * nlev)
        ends = per_level(end, list(si.seqlengths))
        strands = (STAR,) if ignore_strand else (PLUS, MINUS, STAR)
        keys = self._group_keys(ignore_strand)
        pieces = []
        for lv in range(nlev):
            for st in strands:
                idx = np.nonzero(keys == lv * 3 + st)[0]
                gs, gw = _inter.gaps_ranges(self._s()[idx], self._w()[idx], starts[lv], ends[lv])
                if len(gs):
                    pieces.append((lv, st, gs, gw))
        return self._assemble(pieces)

    # ---- set operations -----------------------------------------------------------------------------

    def union(self, y: "GRanges", ignore_strand: bool = False) -> "GRanges":
        return self.concat(y).reduce(drop_empty_ranges=True, ignore_strand=ignore_strand)

    def _binary_groups(self, y: "GRanges", ignore_strand: bool) -> tuple[Seqinfo, dict[int, IntArray], dict[int, IntArray]]:
        si = self._seqinfo.merge(y._seqinfo)
        xk = self._recode_keys(si, ignore_strand)
        yk = y._recode_keys(si, ignore_strand)
        gx = {int(k): np.nonzero(xk == k)[0] for k in np.unique(xk)}
        gy = {int(k): np.nonzero(yk == k)[0] for k in np.unique(yk)}
        return si, gx, gy

    def _recode_keys(self, si: Seqinfo, ignore_strand: bool) -> IntArray:
        st = np.full(len(self), STAR, dtype=np.int64) if ignore_strand else self._strand.astype(np.int64)
        return self._recode(si).astype(np.int64) * 3 + st

    def intersect(self, y: "GRanges", ignore_strand: bool = False) -> "GRanges":
        si, gx, gy = self._binary_groups(y, ignore_strand)
        pieces = []
        for k in sorted(set(gx) & set(gy)):
            s, w = _inter.intersect_ranges(self._s()[gx[k]], self._w()[gx[k]], y._s()[gy[k]], y._w()[gy[k]])
            if len(s):
                pieces.append((k // 3, k % 3, s, w))
        return self._assemble(pieces, seqinfo=si)

    def setdiff(self, y: "GRanges", ignore_strand: bool = False) -> "GRanges":
        si, gx, gy = self._binary_groups(y, ignore_strand)
        pieces = []
        for k in sorted(gx):
            ys = y._s()[gy[k]] if k in gy else np.zeros(0, dtype=np.int64)
            yw = y._w()[gy[k]] if k in gy else np.zeros(0, dtype=np.int64)
            s, w = _inter.setdiff_ranges(self._s()[gx[k]], self._w()[gx[k]], ys, yw)
            if len(s):
                pieces.append((k // 3, k % 3, s, w))
        return self._assemble(pieces, seqinfo=si)

    def _parallel_setup(
        self, y: "GRanges", ignore_strand: bool, allow_recycle_y: bool = False
    ) -> tuple[Seqinfo, "GRanges", npt.NDArray[np.bool_], npt.NDArray[np.int8]]:
        """Merged seqinfo, y aligned to len(self), per-pair compatibility, per-pair result strand."""
        si = self._seqinfo.merge(y._seqinfo)
        if len(y) != len(self):
            if allow_recycle_y and len(y) == 1:
                y = y._take(np.zeros(len(self), dtype=np.int64))
            else:
                raise ValueError("x and y must have the same length")
        same_seq = self._recode(si) == y._recode(si)
        ok = same_seq if ignore_strand else same_seq & _compat(self._strand, y._strand)
        if ignore_strand:
            strand = self._strand.copy()
        else:
            strand = np.where(self._strand == STAR, y._strand, self._strand).astype(np.int8)
        return si, y, ok, strand

    def punion(self, y: "GRanges", fill_gap: bool = False, ignore_strand: bool = False) -> "GRanges":
        si, y, ok, _ = self._parallel_setup(y, ignore_strand)
        if not np.all(ok):
            raise ValueError("x and y elements must have compatible seqnames and strand values")
        r = self._ranges.punion(y._ranges, fill_gap)
        return GRanges._make(self._recode(si), self._strand, r, si)

    def pintersect(
        self,
        y: "GRanges",
        drop_nohit_ranges: bool = False,
        ignore_strand: bool = False,
        resolve_empty: str = "none",
    ) -> "GRanges":
        si, y, ok, strand = self._parallel_setup(y, ignore_strand, allow_recycle_y=True)
        xs, xw = self._s(), self._w()
        s = np.maximum(xs, y._s())
        e = np.minimum(xs + xw - 1, y._s() + y._w() - 1)
        w = e - s + 1
        bad_overlap = ok & (w < 0)
        if bad_overlap.any():
            tmp = IRanges._make(xs[bad_overlap], xw[bad_overlap]).pintersect(
                IRanges._make(y._s()[bad_overlap], y._w()[bad_overlap]), resolve_empty
            )
            s[bad_overlap], w[bad_overlap] = tmp._s(), tmp._w()
        out_s = np.where(ok, s, xs)
        out_w = np.where(ok, w, 0)
        out_strand = np.where(ok, strand, self._strand)
        hit = pl.Series("hit", ok.tolist(), dtype=pl.Boolean)
        base = self._ranges._mcols
        mcols = pl.DataFrame([hit]) if base is None else base.with_columns(hit)
        res = GRanges._make(self._recode(si), out_strand, IRanges._make(out_s, out_w, self._ranges._names, mcols), si)
        if drop_nohit_ranges:
            kept = res._take(np.nonzero(ok)[0])
            return kept.with_mcols(None if base is None else kept.mcols.drop("hit"))
        return res

    def psetdiff(self, y: "GRanges", ignore_strand: bool = False) -> "GRanges":
        si, y, ok, strand = self._parallel_setup(y, ignore_strand)
        xs, xw = self._s(), self._w()
        out_s, out_w = xs.copy(), xw.copy()
        if ok.any():
            sel = np.nonzero(ok)[0]
            r = IRanges._make(xs[sel], xw[sel]).psetdiff(IRanges._make(y._s()[sel], y._w()[sel]))
            out_s[sel], out_w[sel] = r._s(), r._w()
        out_strand = np.where(ok, strand, self._strand)
        return GRanges._make(self._recode(si), out_strand, IRanges._make(out_s, out_w, self._ranges._names), si)

    # ---- overlaps ---------------------------------------------------------------------------------------

    def _overlap_pairs(
        self, subject: "GRanges", maxgap: int, minoverlap: int, type: str, ignore_strand: bool
    ) -> tuple[IntArray, IntArray]:
        _overlaps.check_overlap_args(maxgap, minoverlap, type)
        si = self._seqinfo.merge(subject._seqinfo)
        qc, sc = self._recode(si), subject._recode(si)
        qs_all: list[IntArray] = []
        ss_all: list[IntArray] = []
        for lv in np.intersect1d(qc, sc):
            qi = np.nonzero(qc == lv)[0]
            sj = np.nonzero(sc == lv)[0]
            pq, ps = _overlaps.overlap_pairs(
                self._s()[qi], self._w()[qi], subject._s()[sj], subject._w()[sj], maxgap, minoverlap, type
            )
            q, s = qi[pq], sj[ps]
            if not ignore_strand and len(q):
                keep = _compat(self._strand[q], subject._strand[s])
                q, s = q[keep], s[keep]
            qs_all.append(q)
            ss_all.append(s)
        if not qs_all:
            z = np.zeros(0, dtype=np.int64)
            return z, z
        q, s = np.concatenate(qs_all), np.concatenate(ss_all)
        order = np.lexsort((s, q))
        return q[order], s[order]

    def find_overlaps(
        self,
        subject: "GRanges",
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",
        select: str = "all",
        ignore_strand: bool = False,
    ) -> Hits | npt.NDArray[np.int64]:
        if select not in ("all", "first", "last", "arbitrary"):
            raise ValueError("select must be 'all', 'first', 'last' or 'arbitrary'")
        q, s = self._overlap_pairs(subject, maxgap, minoverlap, type, ignore_strand)
        return _select_hits(q, s, len(self), len(subject), select)

    def count_overlaps(
        self, subject: "GRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any", ignore_strand: bool = False
    ) -> npt.NDArray[np.int64]:
        q, _ = self._overlap_pairs(subject, maxgap, minoverlap, type, ignore_strand)
        return np.bincount(q, minlength=len(self)).astype(np.int64)

    def overlaps_any(
        self, subject: "GRanges", maxgap: int = -1, minoverlap: int = 0, type: str = "any", ignore_strand: bool = False
    ) -> npt.NDArray[np.bool_]:
        return self.count_overlaps(subject, maxgap, minoverlap, type, ignore_strand) > 0

    def subset_by_overlaps(
        self,
        ranges: "GRanges",
        maxgap: int = -1,
        minoverlap: int = 0,
        type: str = "any",
        invert: bool = False,
        ignore_strand: bool = False,
    ) -> "GRanges":
        hit = self.overlaps_any(ranges, maxgap, minoverlap, type, ignore_strand)
        return self[~hit if invert else hit]

    # ---- nearest ---------------------------------------------------------------------------------------------

    def _pair_masks(self, subject: "GRanges", ignore_strand: bool) -> tuple[Any, Any, "GRanges"]:
        si = self._seqinfo.merge(subject._seqinfo)
        qc, sc = self._recode(si), subject._recode(si)

        def mask(i: int) -> npt.NDArray[np.bool_]:
            ok = sc == qc[i]
            return ok if ignore_strand else ok & _compat(self._strand[i], subject._strand)

        def reverse(i: int) -> npt.NDArray[np.bool_]:
            if ignore_strand:
                return np.zeros(len(subject), dtype=bool)
            q = self._strand[i]
            if q == MINUS:
                return np.ones(len(subject), dtype=bool)
            if q == STAR:
                return subject._strand == MINUS
            return np.zeros(len(subject), dtype=bool)

        return mask, reverse, subject

    def nearest(
        self, subject: "GRanges | None" = None, select: str = "arbitrary", ignore_strand: bool = False
    ) -> Hits | npt.NDArray[np.int64]:
        if select not in ("arbitrary", "all"):
            raise ValueError("select must be 'arbitrary' or 'all'")
        me = subject is None
        sj = self if subject is None else subject
        mask, _, _ = self._pair_masks(sj, ignore_strand)
        ties = _nearest.nearest_ties(self._s(), self._w(), sj._s(), sj._w(), mask, me)
        return _ties_result(ties, len(self), len(sj), select)

    def precede(
        self, subject: "GRanges | None" = None, select: str = "first", ignore_strand: bool = False
    ) -> Hits | npt.NDArray[np.int64]:
        if select not in ("first", "all"):
            raise ValueError("select must be 'first' or 'all'")
        return self._directional(subject, "precede", select, ignore_strand)

    def follow(
        self, subject: "GRanges | None" = None, select: str = "last", ignore_strand: bool = False
    ) -> Hits | npt.NDArray[np.int64]:
        if select not in ("last", "all"):
            raise ValueError("select must be 'last' or 'all'")
        return self._directional(subject, "follow", select, ignore_strand)

    def _directional(
        self, subject: "GRanges | None", kind: str, select: str, ignore_strand: bool
    ) -> Hits | npt.NDArray[np.int64]:
        me = subject is None
        sj = self if subject is None else subject
        mask, reverse, _ = self._pair_masks(sj, ignore_strand)
        ties = _nearest.directional_ties(self._s(), self._w(), sj._s(), sj._w(), kind, reverse, mask, me)
        return _ties_result(ties, len(self), len(sj), select)

    def distance(self, other: "GRanges", ignore_strand: bool = False) -> npt.NDArray[np.int64]:
        """Parallel distance; NA_INDEX (-1) where seqnames differ or strands are incompatible."""
        n = max(len(self), len(other)) if len(self) and len(other) else 0
        x = self._take(rep_len(np.arange(len(self)), n)) if n else self[[]]
        y = other._take(rep_len(np.arange(len(other)), n)) if n else other[[]]
        si = x._seqinfo.merge(y._seqinfo)
        ok = x._recode(si) == y._recode(si)
        if not ignore_strand:
            ok &= _compat(x._strand, y._strand)
        d = _nearest.pair_distance(x._s(), x._s() + x._w() - 1, y._s(), y._s() + y._w() - 1)
        return np.where(ok, d, NA_INDEX).astype(np.int64)

    # ---- coverage -------------------------------------------------------------------------------------------------

    def coverage(
        self, shift: Any = 0, width: int | Sequence[int] | Mapping[str, int] | None = None, weight: Any = 1
    ) -> dict[str, Rle]:
        """Per-seqlevel run-length coverage (strand ignored), as {seqlevel: Rle}."""
        si = self._seqinfo
        n = len(self)
        sh = rep_len(as_int64(shift, "shift"), n, "shift") if n else as_int64(shift, "shift")
        if isinstance(weight, str):
            if weight not in self.mcols.columns:
                raise KeyError(f"no mcols column {weight!r}")
            wt = self.mcols[weight].to_numpy()
            if wt.dtype.kind not in "iuf":
                raise TypeError("weight column must be numeric")
        else:
            wt = np.asarray(weight)
            if wt.ndim == 0:
                wt = wt.reshape(1)
            if wt.dtype.kind not in "iuf":
                raise TypeError("weight must be numeric")
        if len(wt) > max(n, 1):
            raise ValueError("weight is longer than x")
        wt = rep_len(wt, n, "weight") if n else wt
        if width is None:
            lens: list[int | None] = list(si.seqlengths)
        elif isinstance(width, Mapping):
            lens = [width[nm] for nm in si.seqnames]
        elif isinstance(width, (int, np.integer)):
            lens = [int(width)] * len(si)
        else:
            lens = [int(v) for v in rep_len(np.asarray(list(width)), len(si), "width")] if len(si) else []
        out: dict[str, Rle] = {}
        for lv, name in enumerate(si.seqnames):
            idx = np.nonzero(self._seq == lv)[0]
            length = lens[lv]
            if length is None:
                ends = self._s()[idx] + self._w()[idx] - 1 + (sh[idx] if n else 0)
                length = int(max(ends.max(), 0)) if len(idx) else 0
            out[name] = _coverage.coverage_rle(
                self._s()[idx], self._w()[idx], length, sh[idx] if n else sh, wt[idx] if n else wt
            )
        return out

    # ---- ordering ---------------------------------------------------------------------------------------------------

    def _keys(self, ignore_strand: bool = False) -> list[IntArray]:
        keys = [self._seq.astype(np.int64)]
        if not ignore_strand:
            keys.append(self._strand.astype(np.int64))
        return keys + [self._s(), self._w()]

    def order(self, decreasing: bool = False, ignore_strand: bool = False) -> IntArray:
        return _ordering.order_by(self._keys(ignore_strand), decreasing)

    def sort(self, decreasing: bool = False, ignore_strand: bool = False) -> "GRanges":
        return self._take(self.order(decreasing, ignore_strand))

    def rank(self, ties_method: str = "average", ignore_strand: bool = False) -> npt.NDArray[np.generic]:
        return _ordering.rank_by(self._keys(ignore_strand), ties_method)

    def duplicated(self, from_last: bool = False) -> npt.NDArray[np.bool_]:
        return _ordering.duplicated_by(self._keys(), from_last)

    def unique(self, from_last: bool = False) -> "GRanges":
        return self[~self.duplicated(from_last)]

    # ---- misc -----------------------------------------------------------------------------------------------------------

    def __repr__(self) -> str:
        rows = [
            f"{sn}:{s}-{s + w - 1}:{st}"
            for sn, s, w, st in zip(
                self.seqname_strings()[:5], self._s()[:5].tolist(), self._w()[:5].tolist(), self.strand_strings()[:5]
            )
        ]
        return f"GRanges(n={len(self)}: {', '.join(rows)}{', ...' if len(self) > 5 else ''})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, GRanges):
            return NotImplemented
        return (
            np.array_equal(self._seq, other._seq)
            and np.array_equal(self._strand, other._strand)
            and self._ranges == other._ranges
            and self._seqinfo == other._seqinfo
        )

    __hash__ = None  # type: ignore[assignment]


def _series_values(v: Any) -> Any:
    if isinstance(v, (pl.Series, list)):
        return v
    return v.tolist() if isinstance(v, np.ndarray) else v
