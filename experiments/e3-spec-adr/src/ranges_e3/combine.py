"""``concat``: append IRanges (and, later, GRanges) in order."""

from __future__ import annotations

from ranges_e3.iranges import IRanges, concat_iranges


def concat(*objs: IRanges) -> IRanges:
    """Append objects of one type in argument order (CORE.11, GR.5)."""
    if not objs:
        raise ValueError("concat needs at least one object")
    if all(isinstance(o, IRanges) for o in objs):
        return concat_iranges(objs)
    raise TypeError("concat needs IRanges objects of one type")
