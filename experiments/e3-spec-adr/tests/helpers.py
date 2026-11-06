"""Small constructors/readers so tests can state ranges as ``[s,e]`` pairs."""

from __future__ import annotations

from typing import Any

import numpy as np

from ranges_e3 import IRanges


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
