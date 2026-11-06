"""``concat``: append objects of one type in order."""

from __future__ import annotations

from typing import TypeVar

from ranges_e3.granges import GRanges, concat_granges
from ranges_e3.iranges import IRanges, concat_iranges

T = TypeVar("T", IRanges, GRanges)


def concat(*objs: T) -> T:
    """Append objects in argument order (CORE.11, GR.5); types must agree."""
    if not objs:
        raise ValueError("concat needs at least one object")
    if all(isinstance(o, IRanges) for o in objs):
        return concat_iranges(objs)  # type: ignore[return-value, arg-type]
    if all(isinstance(o, GRanges) for o in objs):
        return concat_granges(objs)  # type: ignore[return-value, arg-type]
    raise TypeError("concat needs objects of one type (IRanges or GRanges)")
