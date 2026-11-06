"""Seqinfo: per-sequence metadata (name, length, circularity, genome)."""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence


class Seqinfo:
    """Immutable table of sequence levels, in a meaningful order.

    ``seqlengths``, ``is_circular`` and ``genome`` use ``None`` for R's ``NA``.
    """

    __slots__ = ("_names", "_lengths", "_circular", "_genome", "_index")

    def __init__(
        self,
        seqnames: Sequence[str],
        seqlengths: Sequence[int | None] | None = None,
        is_circular: Sequence[bool | None] | None = None,
        genome: str | Sequence[str | None] | None = None,
    ) -> None:
        names = tuple(str(s) for s in seqnames)
        if len(set(names)) != len(names):
            raise ValueError("'seqnames' must not contain duplicates")
        n = len(names)
        self._names = names
        self._lengths = _fill(seqlengths, n, "seqlengths", int)
        self._circular = _fill(is_circular, n, "is_circular", bool)
        if isinstance(genome, str):
            genome = [genome] * n
        self._genome = _fill(genome, n, "genome", str)
        for length in self._lengths:
            if length is not None and length < 0:
                raise ValueError("seqlengths must be non-negative")
        self._index = {s: i for i, s in enumerate(names)}

    # ----- accessors ------------------------------------------------------
    @property
    def seqnames(self) -> tuple[str, ...]:
        return self._names

    seqlevels = seqnames

    @property
    def seqlengths(self) -> dict[str, int | None]:
        return dict(zip(self._names, self._lengths))

    @property
    def is_circular(self) -> dict[str, bool | None]:
        return dict(zip(self._names, self._circular))

    @property
    def genome(self) -> dict[str, str | None]:
        return dict(zip(self._names, self._genome))

    def __len__(self) -> int:
        return len(self._names)

    def index(self, name: str) -> int:
        """Position of a sequence level (``KeyError`` if absent)."""
        return self._index[name]

    def __contains__(self, name: object) -> bool:
        return name in self._index

    # ----- combination ----------------------------------------------------
    def merge(self, other: Seqinfo) -> Seqinfo:
        """R ``merge``: union of the levels (this object's order first, then new
        levels of ``other``). ``None`` fields are filled from the other side;
        conflicting non-``None`` lengths, circularity flags or genomes raise
        ``ValueError``. Warns when either side has levels the other lacks."""
        only_self = [s for s in self._names if s not in other._index]
        only_other = [s for s in other._names if s not in self._index]
        if only_self and only_other:
            warnings.warn(
                "Each of the 2 combined objects has sequence levels not in the other: "
                f"in 'x': {', '.join(only_self)}; in 'y': {', '.join(only_other)}",
                stacklevel=3,
            )
        names = list(self._names) + only_other
        lengths = list(self._lengths) + [None] * len(only_other)
        circ = list(self._circular) + [None] * len(only_other)
        genome = list(self._genome) + [None] * len(only_other)
        for j, name in enumerate(other._names):
            i = names.index(name)
            for col, mine, theirs, what in (
                (lengths, lengths[i], other._lengths[j], "seqlengths"),
                (circ, circ[i], other._circular[j], "circularity flags"),
                (genome, genome[i], other._genome[j], "genomes"),
            ):
                if mine is None:
                    col[i] = theirs
                elif theirs is not None and mine != theirs:
                    raise ValueError(
                        f"sequence {name} has incompatible {what}: in 'x': {mine}; in 'y': {theirs}"
                    )
        return Seqinfo(names, lengths, circ, genome)

    # ----- comparison / display --------------------------------------------
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Seqinfo):
            return NotImplemented
        return (
            self._names == other._names
            and self._lengths == other._lengths
            and self._circular == other._circular
            and self._genome == other._genome
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        rows = ", ".join(
            f"{n}(len={l}, circ={c}, genome={g})"
            for n, l, c, g in zip(self._names, self._lengths, self._circular, self._genome)
        )
        return f"Seqinfo([{rows}])"

    @classmethod
    def from_mapping(cls, seqlengths: Mapping[str, int | None]) -> Seqinfo:
        """Build from ``{seqname: length}`` (R's ``seqlengths=`` named vector)."""
        return cls(list(seqlengths), list(seqlengths.values()))


def _fill(values: Sequence | None, n: int, name: str, cast: type) -> tuple:
    if values is None:
        return (None,) * n
    vals = list(values)
    if len(vals) != n:
        raise ValueError(f"'{name}' must have one entry per sequence ({n}), got {len(vals)}")
    return tuple(None if v is None else cast(v) for v in vals)
