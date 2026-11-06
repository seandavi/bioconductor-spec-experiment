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

from . import _order, _validate as V
from .iranges import IRanges
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
        seq = np.concatenate([_recode(o._seq, o._seqinfo, si) for o in objs]) if len(objs) else None
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
