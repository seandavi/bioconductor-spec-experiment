"""Seqinfo: per-sequence name, length, circularity, genome (None = R NA)."""

import warnings
from collections.abc import Sequence
from typing import Any


def _field(value: Any, n: int, name: str) -> tuple[Any, ...]:
    if value is None:
        return (None,) * n
    if isinstance(value, (str, bool, int)):
        return (value,) * n
    out = tuple(value)
    if len(out) != n:
        raise ValueError(f"{name} must have length {n}, got {len(out)}")
    return out


class Seqinfo:
    """Immutable table of sequences; order defines the seqlevels order."""

    __slots__ = ("_seqnames", "_seqlengths", "_is_circular", "_genome")

    def __init__(
        self,
        seqnames: Sequence[str],
        seqlengths: Any = None,
        is_circular: Any = None,
        genome: Any = None,
    ) -> None:
        names = tuple(seqnames)
        if len(set(names)) != len(names):
            raise ValueError("seqnames must be unique")
        n = len(names)
        lengths = _field(seqlengths, n, "seqlengths")
        if any(v is not None and (not isinstance(v, int) or v < 0) for v in lengths):
            raise ValueError("seqlengths must be non-negative integers or None")
        circ = _field(is_circular, n, "is_circular")
        if any(v is not None and not isinstance(v, bool) for v in circ):
            raise ValueError("is_circular must be bool or None")
        self._seqnames = names
        self._seqlengths = tuple(None if v is None else int(v) for v in lengths)
        self._is_circular = circ
        self._genome = _field(genome, n, "genome")

    @property
    def seqnames(self) -> tuple[str, ...]:
        return self._seqnames

    @property
    def seqlengths(self) -> tuple[int | None, ...]:
        return self._seqlengths

    @property
    def is_circular(self) -> tuple[bool | None, ...]:
        return self._is_circular

    @property
    def genome(self) -> tuple[str | None, ...]:
        return self._genome

    def __len__(self) -> int:
        return len(self._seqnames)

    def _key(self) -> tuple[Any, ...]:
        return (self._seqnames, self._seqlengths, self._is_circular, self._genome)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Seqinfo) and self._key() == other._key()

    def __hash__(self) -> int:
        return hash(self._key())

    def __repr__(self) -> str:
        rows = [
            f"  {nm}: length={ln}, circular={ci}, genome={ge}"
            for nm, ln, ci, ge in zip(self._seqnames, self._seqlengths, self._is_circular, self._genome)
        ]
        return f"Seqinfo of {len(self)} sequences\n" + "\n".join(rows)

    def index(self, name: str) -> int:
        try:
            return self._seqnames.index(name)
        except ValueError:
            raise KeyError(f"unknown seqlevel {name!r}") from None

    def merge(self, other: "Seqinfo") -> "Seqinfo":
        """Union of sequences: this one's order, then new ones from `other`.

        Per field, None yields to a value; two different values raise ValueError.
        """
        if self._seqnames and other._seqnames and not set(self._seqnames) & set(other._seqnames):
            warnings.warn("the two Seqinfo objects have no sequence names in common", stacklevel=2)
        names = list(self._seqnames) + [n for n in other._seqnames if n not in set(self._seqnames)]
        fields = {"seqlengths": [], "is_circular": [], "genome": []}  # type: dict[str, list[Any]]
        for nm in names:
            for attr, out in fields.items():
                vals = []
                for si in (self, other):
                    if nm in si._seqnames:
                        vals.append(getattr(si, attr)[si._seqnames.index(nm)])
                known = {v for v in vals if v is not None}
                if len(known) > 1:
                    raise ValueError(f"conflicting {attr} for sequence {nm!r}: {sorted(known, key=str)}")
                out.append(known.pop() if known else None)
        return Seqinfo(names, fields["seqlengths"], fields["is_circular"], fields["genome"])
