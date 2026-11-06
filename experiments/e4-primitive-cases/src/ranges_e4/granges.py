"""GRanges: IRanges + seqnames + strand + seqinfo."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Self

import numpy as np
import polars as pl

from . import _mcols, intra
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
