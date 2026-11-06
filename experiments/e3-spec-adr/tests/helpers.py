"""Small constructors/readers so tests can state ranges as ``[s,e]`` pairs."""

from __future__ import annotations

from typing import Any

import numpy as np

from ranges_e3 import GRanges, IRanges


def ir(*pairs: tuple[int, int], **kwargs: Any) -> IRanges:
    """``ir((1,3),(5,4))`` -> ``{[1,3],[5,4]}``."""
    if not pairs:
        return IRanges(start=[], width=[], **kwargs)
    return IRanges(start=[p[0] for p in pairs], end=[p[1] for p in pairs], **kwargs)


def se(x: IRanges) -> list[tuple[int, int]]:
    """Ranges as ``[(start, end), ...]``."""
    return [(int(s), int(e)) for s, e in zip(x.start, x.end)]


def na(arr: np.ma.MaskedArray[Any, Any]) -> list[int | None]:
    """Masked int array as a list with ``None`` for NA."""
    return [None if m else int(v) for v, m in zip(arr.data, np.ma.getmaskarray(arr))]


def gr(*rows: tuple[str, int, int, str], **kwargs: Any) -> GRanges:
    """``gr(("chr1",1,3,"+"), ...)`` -> GRanges (seqname, start, end, strand)."""
    if not rows:
        return GRanges([], IRanges(start=[], width=[]), **kwargs)
    return GRanges(
        [r[0] for r in rows],
        IRanges(start=[r[1] for r in rows], end=[r[2] for r in rows]),
        strand=[r[3] for r in rows],
        **kwargs,
    )


def gs(x: GRanges) -> list[tuple[str, int, int, str]]:
    """GRanges as ``[(seqname, start, end, strand), ...]``."""
    return [
        (str(a), int(b), int(c), str(d)) for a, b, c, d in zip(x.seqnames, x.start, x.end, x.strand)
    ]
