"""Seqinfo: per-sequence name, length, circularity, genome (None = NA)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any


def _broadcast(values: Any, n: int, what: str) -> tuple[Any, ...]:
    if values is None:
        return (None,) * n
    if isinstance(values, (str, int, bool)):
        values = [values]
    values = tuple(values)
    if len(values) == 1 and n != 1:
        return values * n
    if len(values) != n:
        raise ValueError(f"'{what}' has length {len(values)}, expected 1 or {n}")
    return values


@dataclass(frozen=True)
class Seqinfo:
    seqnames: tuple[str, ...]
    _seqlengths: tuple[int | None, ...]
    _is_circular: tuple[bool | None, ...]
    _genome: tuple[str | None, ...]

    def __init__(
        self,
        seqnames: Sequence[str],
        seqlengths: Any = None,
        is_circular: Any = None,
        genome: Any = None,
    ) -> None:
        names = tuple(seqnames)
        if any(s is None or s == "" for s in names):
            raise ValueError("seqnames cannot contain NA or empty strings")
        if len(set(names)) != len(names):
            raise ValueError("seqnames must be unique")
        lengths = _broadcast(seqlengths, len(names), "seqlengths")
        if any(v is not None and v < 0 for v in lengths):
            raise ValueError("seqlengths must be non-negative")
        object.__setattr__(self, "seqnames", names)
        object.__setattr__(self, "_seqlengths", tuple(None if v is None else int(v) for v in lengths))
        object.__setattr__(
            self, "_is_circular", tuple(None if v is None else bool(v) for v in _broadcast(is_circular, len(names), "is_circular"))
        )
        object.__setattr__(self, "_genome", _broadcast(genome, len(names), "genome"))

    def __len__(self) -> int:
        return len(self.seqnames)

    # named-vector accessors (dict: seqname -> value) ----------------------
    @property
    def seqlevels(self) -> list[str]:
        return list(self.seqnames)

    @property
    def seqlengths(self) -> dict[str, int | None]:
        return dict(zip(self.seqnames, self._seqlengths, strict=True))

    @property
    def is_circular(self) -> dict[str, bool | None]:
        return dict(zip(self.seqnames, self._is_circular, strict=True))

    @property
    def genome(self) -> dict[str, str | None]:
        return dict(zip(self.seqnames, self._genome, strict=True))

    def index(self, name: str) -> int:
        return self.seqnames.index(name)

    def merge(self, other: Seqinfo) -> Seqinfo:
        """Union by name (self's levels first). Non-NA values fill NA; conflicts raise."""
        names = list(self.seqnames) + [n for n in other.seqnames if n not in self.seqnames]
        lengths, circ, genome = [], [], []
        for n in names:
            vals = []
            for si in (self, other):
                if n in si.seqnames:
                    i = si.index(n)
                    vals.append((si._seqlengths[i], si._is_circular[i], si._genome[i]))
            for what, k, out in (("seqlengths", 0, lengths), ("circularity flags", 1, circ), ("genomes", 2, genome)):
                known = {v[k] for v in vals if v[k] is not None}
                if len(known) > 1:
                    raise ValueError(f"sequence {n} has incompatible {what}")
                out.append(next(iter(known)) if known else None)
        return Seqinfo(names, lengths, circ, genome)
