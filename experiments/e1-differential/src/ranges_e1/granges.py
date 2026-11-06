"""GRanges: IRanges plus seqnames, strand and seqinfo."""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _meta
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
