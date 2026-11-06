"""Seqinfo: per-sequence name, length, circularity and genome, with R's merge rules."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any


def _broadcast(value: Any, n: int, what: str) -> tuple[Any, ...]:
    if value is None:
        return (None,) * n
    if isinstance(value, (str, bool, int)) or not isinstance(value, Sequence):
        value = [value]
    seq = tuple(value)
    if len(seq) == 1 and n != 1:
        seq = seq * n
    if len(seq) != n:
        raise ValueError(f"{what} must have length {n} (the number of sequences), got {len(seq)}")
    return seq


@dataclass(frozen=True)
class Seqinfo:
    """Immutable table of sequences. `None` is NA in any field."""

    seqnames: tuple[str, ...]
    seqlengths: tuple[int | None, ...]
    is_circular: tuple[bool | None, ...]
    genome: tuple[str | None, ...]

    def __init__(
        self,
        seqnames: Sequence[str] = (),
        seqlengths: Sequence[int | None] | int | None = None,
        is_circular: Sequence[bool | None] | bool | None = None,
        genome: Sequence[str | None] | str | None = None,
    ) -> None:
        names = tuple(seqnames)
        if any(n is None or n == "" for n in names):
            raise ValueError("seqnames cannot contain NAs or empty strings")
        if len(set(names)) != len(names):
            raise ValueError("seqnames cannot contain duplicated sequence names")
        lens = _broadcast(seqlengths, len(names), "seqlengths")
        if any(v is not None and (int(v) != v or v < 0) for v in lens):
            raise ValueError("seqlengths must be non-negative integers")
        object.__setattr__(self, "seqnames", names)
        object.__setattr__(self, "seqlengths", tuple(None if v is None else int(v) for v in lens))
        object.__setattr__(self, "is_circular", _broadcast(is_circular, len(names), "is_circular"))
        object.__setattr__(self, "genome", _broadcast(genome, len(names), "genome"))

    def __len__(self) -> int:
        return len(self.seqnames)

    def index_of(self, name: str) -> int:
        try:
            return self.seqnames.index(name)
        except ValueError:
            raise KeyError(f"unknown sequence {name!r}") from None

    def subset(self, names: Sequence[str]) -> "Seqinfo":
        """Select/reorder levels by name; unknown names add a level with NA metadata."""
        rows = []
        for nm in names:
            if nm in self.seqnames:
                i = self.seqnames.index(nm)
                rows.append((nm, self.seqlengths[i], self.is_circular[i], self.genome[i]))
            else:
                rows.append((nm, None, None, None))
        return Seqinfo([r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows], [r[3] for r in rows])

    def merge(self, other: "Seqinfo") -> "Seqinfo":
        """Union of levels (self's order first); NA fields filled from the other side; conflicts error."""
        names = list(self.seqnames)
        lens = list(self.seqlengths)
        circ = list(self.is_circular)
        genome = list(self.genome)
        for j, nm in enumerate(other.seqnames):
            if nm not in self.seqnames:
                names.append(nm)
                lens.append(other.seqlengths[j])
                circ.append(other.is_circular[j])
                genome.append(other.genome[j])
                continue
            i = self.seqnames.index(nm)
            fields: list[tuple[str, Any, Any, list[Any]]] = [
                ("seqlengths", self.seqlengths[i], other.seqlengths[j], lens),
                ("circularity flags", self.is_circular[i], other.is_circular[j], circ),
                ("genomes", self.genome[i], other.genome[j], genome),
            ]
            for label, mine, theirs, target in fields:
                if mine is not None and theirs is not None and mine != theirs:
                    raise ValueError(f"sequence {nm} has incompatible {label}: {mine} vs {theirs}")
                if mine is None:
                    target[i] = theirs
        return Seqinfo(names, lens, circ, genome)

    def with_field(self, field: str, value: Any) -> "Seqinfo":
        """Return a copy with seqlengths / is_circular / genome replaced (scalar recycles)."""
        n = len(self)
        if field == "seqlengths":
            return Seqinfo(self.seqnames, _broadcast(value, n, "seqlengths"), self.is_circular, self.genome)
        if field == "is_circular":
            return Seqinfo(self.seqnames, self.seqlengths, _broadcast(value, n, "is_circular"), self.genome)
        if field == "genome":
            return Seqinfo(self.seqnames, self.seqlengths, self.is_circular, _broadcast(value, n, "genome"))
        raise ValueError(f"unknown field {field}")

    def rename_levels(self, value: Sequence[str]) -> "Seqinfo":
        """R's `seqlevels<-` on a Seqinfo: positional rename if nothing in common, else map by name.

        Returns the new Seqinfo; use `level_map` for the old->new index mapping on data.
        """
        new, _ = self.relevel(value)
        return new

    def relevel(self, value: Sequence[str]) -> tuple["Seqinfo", list[int | None]]:
        """Return (new Seqinfo, old_index -> new_index or None when the level is dropped)."""
        new_names = list(value)
        if len(set(new_names)) != len(new_names):
            raise ValueError("new seqlevels cannot contain duplicates")
        common = [v for v in new_names if v in self.seqnames]
        if not common and len(new_names) == len(self.seqnames):
            mapping: list[int | None] = list(range(len(self.seqnames)))
            return (
                Seqinfo(new_names, self.seqlengths, self.is_circular, self.genome),
                mapping,
            )
        new = self.subset(new_names)
        mapping = [new_names.index(nm) if nm in new_names else None for nm in self.seqnames]
        return new, mapping
