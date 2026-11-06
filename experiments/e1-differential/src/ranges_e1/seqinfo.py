"""Seqinfo: per-sequence name, length, circularity and genome."""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np


def _per_sequence(v: Any, n: int, what: str, check: Any, *, recycle: bool = False) -> tuple[Any, ...]:
    """Normalise an optional per-sequence argument to a tuple of length n (None = NA).

    R recycles a single `genome` only; seqlengths / isCircular must have exactly n values.
    """
    if v is None:
        return (None,) * n
    if isinstance(v, (str, bool, int, float)):
        v = [v]
    if np.ma.isMaskedArray(v):
        v = [None if m else x for x, m in zip(np.ma.getdata(v).tolist(), np.ma.getmaskarray(v))]
    elif isinstance(v, np.ndarray):
        v = v.tolist()
    vals = [None if (x is None or (isinstance(x, float) and x != x)) else check(x) for x in list(v)]
    if recycle and len(vals) == 1 and n != 1:
        vals = vals * n
    if len(vals) != n:
        raise ValueError(f"the length of the supplied '{what}' vector must be equal to the number of sequences")
    return tuple(vals)


def _length(x: Any) -> int:
    if isinstance(x, bool):
        raise ValueError("'seqlengths' must be integers")
    x = int(x)  # R truncates doubles (as.integer)
    if x < 0:
        raise ValueError("the supplied 'seqlengths' vector contains negative values")
    return int(x)


@dataclass(frozen=True)
class Seqinfo:
    """Immutable sequence information. `None` entries are NA."""

    seqnames: tuple[str, ...]
    seqlengths: tuple[int | None, ...]
    is_circular: tuple[bool | None, ...]
    genome: tuple[str | None, ...]

    def __init__(
        self,
        seqnames: Iterable[str] = (),
        seqlengths: Any = None,
        is_circular: Any = None,
        genome: Any = None,
    ) -> None:
        names = tuple(seqnames) if not isinstance(seqnames, str) else (seqnames,)
        if any(s is None or s == "" for s in names):
            raise ValueError("supplied 'seqlevels' cannot contain NAs or empty strings (\"\")")
        if len(set(names)) != len(names):
            raise ValueError("supplied 'seqlevels' cannot contain duplicated sequence names")
        n = len(names)
        object.__setattr__(self, "seqnames", names)
        object.__setattr__(self, "seqlengths", _per_sequence(seqlengths, n, "seqlengths", _length))
        object.__setattr__(self, "is_circular", _per_sequence(is_circular, n, "isCircular", bool))
        object.__setattr__(self, "genome", _per_sequence(genome, n, "genome", str, recycle=True))

    def __len__(self) -> int:
        return len(self.seqnames)

    def index(self, name: str) -> int:
        return self.seqnames.index(name)

    def merge(self, other: Seqinfo) -> Seqinfo:
        """Union of sequences (self's order first); NA fills from the other side, conflicts error.

        Warns when both sides hold sequences the other lacks (R does).
        """
        mine, theirs = set(self.seqnames), set(other.seqnames)
        if not (mine & theirs) and (mine or theirs):
            warnings.warn("The 2 combined objects have no sequence levels in common.", UserWarning, stacklevel=3)
        elif mine - theirs and theirs - mine:
            warnings.warn(
                "Each of the 2 combined objects has sequence levels not in the other.", UserWarning, stacklevel=3
            )
        names = list(self.seqnames) + [s for s in other.seqnames if s not in mine]
        lens: list[int | None] = []
        circ: list[bool | None] = []
        gen: list[str | None] = []

        def pick(what: str, name: str, x: Any, y: Any) -> Any:
            if x is not None and y is not None and x != y:
                raise ValueError(f"sequence {name} has incompatible {what}:\n  - in 'x': {x}\n  - in 'y': {y}")
            return x if x is not None else y

        for s in names:
            i = self.index(s) if s in mine else None
            j = other.index(s) if s in theirs else None
            lens.append(pick("seqlengths", s, None if i is None else self.seqlengths[i], None if j is None else other.seqlengths[j]))
            circ.append(pick("is_circular", s, None if i is None else self.is_circular[i], None if j is None else other.is_circular[j]))
            gen.append(pick("genome", s, None if i is None else self.genome[i], None if j is None else other.genome[j]))
        return Seqinfo(names, lens, circ, gen)

    @classmethod
    def merge_all(cls, parts: Sequence[Seqinfo]) -> Seqinfo:
        out = parts[0]
        for p in parts[1:]:
            out = out.merge(p)
        return out
