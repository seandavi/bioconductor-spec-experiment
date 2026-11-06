"""``Seqinfo``: ordered sequence names with length / circularity / genome (spec 90-seqinfo.md)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TypeVar

T = TypeVar("T")


def _column(values: Sequence[T | None] | None, n: int, what: str) -> tuple[T | None, ...]:
    if values is None:
        return (None,) * n
    out = tuple(values)
    if len(out) != n:
        raise ValueError(f"{what} has length {len(out)}, expected {n}")
    return out


class Seqinfo:
    """Immutable (SEQI.1). Attributes are ``None`` when absent."""

    __slots__ = ("_names", "_lengths", "_circular", "_genome", "_index")

    def __init__(
        self,
        seqnames: Sequence[str],
        seqlengths: Sequence[int | None] | None = None,
        is_circular: Sequence[bool | None] | None = None,
        genome: Sequence[str | None] | None = None,
    ) -> None:
        self._names = tuple(seqnames)
        n = len(self._names)
        if len(set(self._names)) != n:
            raise ValueError("seqnames must be distinct")
        self._lengths = _column(seqlengths, n, "seqlengths")
        if any(v is not None and v < 0 for v in self._lengths):
            raise ValueError("seqlengths must be non-negative")
        self._circular = _column(is_circular, n, "is_circular")
        self._genome = _column(genome, n, "genome")
        self._index = {name: i for i, name in enumerate(self._names)}

    @property
    def seqnames(self) -> list[str]:
        return list(self._names)

    @property
    def seqlengths(self) -> list[int | None]:
        return list(self._lengths)

    @property
    def is_circular(self) -> list[bool | None]:
        return list(self._circular)

    @property
    def genome(self) -> list[str | None]:
        return list(self._genome)

    def __len__(self) -> int:
        return len(self._names)

    def index(self, name: str) -> int:
        """Position of ``name`` among the levels (KeyError if absent)."""
        return self._index[name]

    def __contains__(self, name: object) -> bool:
        return name in self._index

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Seqinfo):
            return NotImplemented
        return (self._names, self._lengths, self._circular, self._genome) == (
            other._names, other._lengths, other._circular, other._genome,
        )

    __hash__ = None  # type: ignore[assignment]

    def merge(self, other: Seqinfo) -> Seqinfo:
        """SEQI.2: union of levels (self's first); conflicting present values are refused."""
        names = list(self._names) + [n for n in other._names if n not in self._index]
        merged: dict[str, list[object]] = {}
        for name in names:
            row: list[object] = []
            for attr in ("_lengths", "_circular", "_genome"):
                a = getattr(self, attr)[self._index[name]] if name in self._index else None
                b = getattr(other, attr)[other._index[name]] if name in other._index else None
                if a is not None and b is not None and a != b:
                    raise ValueError(f"sequence {name!r} has incompatible {attr.strip('_')}: {a!r} vs {b!r}")
                row.append(a if a is not None else b)
            merged[name] = row
        return Seqinfo(
            names,
            [merged[n][0] for n in names],  # type: ignore[misc]
            [merged[n][1] for n in names],  # type: ignore[misc]
            [merged[n][2] for n in names],  # type: ignore[misc]
        )

    def __repr__(self) -> str:
        return f"Seqinfo({list(self._names)}, lengths={list(self._lengths)})"
