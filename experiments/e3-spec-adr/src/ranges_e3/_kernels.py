"""Range math on ``(start, width)`` int64 arrays. No classes (deep module, clear seam).

Every function is pure: arrays in, arrays out. Element indices are 0-based.
"""

from __future__ import annotations

import numpy as np

from ranges_e3._types import BoolArray, IntArray

# --- intra-range (spec 20-intra.md) ---

FIX_START, FIX_END, FIX_CENTER = 0, 1, 2


def shift(start: IntArray, width: IntArray, by: IntArray) -> tuple[IntArray, IntArray]:
    """INTRA.1."""
    return start + by, width


def narrow(
    start: IntArray,
    width: IntArray,
    rel_start: IntArray | None,
    rel_end: IntArray | None,
    rel_width: IntArray | None,
) -> tuple[IntArray, IntArray]:
    """INTRA.2. Relative arguments may be negative; refuses results outside the range."""
    given = (rel_start is not None, rel_end is not None, rel_width is not None)
    if sum(given) == 0 or sum(given) == 3 or given == (False, False, True):
        raise ValueError("narrow needs start and/or end, plus optionally width (not width alone)")
    end = start + width - 1
    a = None if rel_start is None else np.where(rel_start < 0, width + rel_start + 1, rel_start)
    b = None if rel_end is None else np.where(rel_end < 0, width + rel_end + 1, rel_end)
    if a is not None and b is not None:
        new_s, new_e = start + a - 1, start + b - 1
    elif a is not None and rel_width is not None:
        new_s = start + a - 1
        new_e = new_s + rel_width - 1
    elif b is not None and rel_width is not None:
        new_e = start + b - 1
        new_s = new_e - rel_width + 1
    elif a is not None:
        new_s, new_e = start + a - 1, end
    else:
        assert b is not None
        new_s, new_e = start, start + b - 1
    if ((new_s < start) | (new_e > end) | (new_e < new_s - 1)).any():
        raise ValueError("narrow result must lie inside each range")
    return new_s, new_e - new_s + 1


def resize(
    start: IntArray, width: IntArray, new_width: IntArray, fix: IntArray
) -> tuple[IntArray, IntArray]:
    """INTRA.3. ``fix`` holds FIX_* codes per range."""
    if (new_width < 0).any():
        raise ValueError("resize width must be non-negative")
    end = start + width - 1
    new_start = np.where(
        fix == FIX_START,
        start,
        np.where(fix == FIX_END, end - new_width + 1, start + (width - new_width) // 2),
    )
    return new_start, new_width


def flank(
    start: IntArray, width: IntArray, k: IntArray, at_start: BoolArray, both: bool
) -> tuple[IntArray, IntArray]:
    """INTRA.4. ``at_start`` selects the start side per range."""
    end = start + width - 1
    if both:
        a = np.abs(k)
        return np.where(at_start, start - a, end - a + 1), 2 * a
    outward = k >= 0
    new_start = np.where(
        at_start,
        np.where(outward, start - k, start),
        np.where(outward, end + 1, end + k + 1),
    )
    return new_start, np.abs(k)


def promoters(
    start: IntArray, width: IntArray, upstream: IntArray, downstream: IntArray, mirror: BoolArray
) -> tuple[IntArray, IntArray]:
    """INTRA.5. ``mirror`` ranges (``-`` strand) get ``[e-down+1, e+up]``."""
    if (upstream < 0).any() or (downstream < 0).any():
        raise ValueError("promoters upstream/downstream must be non-negative")
    end = start + width - 1
    return np.where(mirror, end - downstream + 1, start - upstream), upstream + downstream
