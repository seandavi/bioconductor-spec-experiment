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
    concat_names,
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
        return self._mcols

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
