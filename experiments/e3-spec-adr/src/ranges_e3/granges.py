"""``GRanges``: IRanges + seqnames + strand + seqinfo (spec 95-granges.md)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Self

import numpy as np
import numpy.typing as npt
import polars as pl

from ranges_e3 import _kernels, _mcols, _subscript
from ranges_e3._types import BoolArray, IntArray
from ranges_e3.iranges import IRanges, per_range
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
