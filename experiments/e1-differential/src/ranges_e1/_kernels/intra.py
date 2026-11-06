"""Intra-range ops on (start, width) arrays. All int64; NA never reaches here."""

from __future__ import annotations

import numpy as np

NA = np.iinfo(np.int64).min  # internal sentinel for "argument not supplied"


def shift(start: np.ndarray, width: np.ndarray, by: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return start + by, width


def narrow(
    start: np.ndarray,
    width: np.ndarray,
    sew_start: np.ndarray,
    sew_end: np.ndarray,
    sew_width: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """R's `solveUserSEW` (allow.nonnarrowing=FALSE) applied to ranges.

    `sew_*` are per-range, relative to each range; NA is the `NA` sentinel.
    """
    ref = width
    s, e, w = sew_start.copy(), sew_end.copy(), sew_width.copy()
    has_w = w != NA
    if (has_w & ((s != NA) == (e != NA))).any():
        raise ValueError(
            "either the supplied start or the supplied end (but not both) "
            "must be NA when the supplied width is not NA"
        )
    neg = (s != NA) & (s < 0)
    s = np.where(neg, s + ref + 1, s)
    neg = (e != NA) & (e < 0)
    e = np.where(neg, e + ref + 1, e)
    if (has_w & (w < 0)).any():
        raise ValueError("the supplied width is negative")
    # fill the missing one
    from_w_end = has_w & (s == NA)
    from_w_start = has_w & (e == NA)
    s = np.where(from_w_end, e - w + 1, s)
    e = np.where(from_w_start, s + w - 1, e)
    s = np.where(s == NA, 1, s)
    e = np.where(e == NA, ref, e)
    if (s < 1).any() or (s > ref + 1).any():
        raise ValueError("the supplied start is out of the range's bounds")
    if (e > ref).any() or (e < 0).any():
        raise ValueError("the supplied end is out of the range's bounds")
    if (e - s + 1 < 0).any():
        raise ValueError("the supplied start/end give a negative width")
    return start + s - 1, e - s + 1


FIX_START, FIX_END, FIX_CENTER = 0, 1, 2


def resize(
    start: np.ndarray, width: np.ndarray, new_width: np.ndarray, fix: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """`fix`: per-range code (FIX_START / FIX_END / FIX_CENTER)."""
    end = start + width - 1
    new_start = np.where(
        fix == FIX_START,
        start,
        np.where(fix == FIX_END, end - new_width + 1, start + (width - new_width) // 2),
    )
    return new_start, new_width


def flank(
    start: np.ndarray, width: np.ndarray, by: np.ndarray, at_start: np.ndarray, both: bool
) -> tuple[np.ndarray, np.ndarray]:
    """Flanking regions; `by < 0` flanks inside the range. `at_start` is per-range."""
    end = start + width - 1
    a = np.abs(by)
    if both:
        new_start = np.where(at_start, start - a, end - a + 1)
        return new_start, 2 * a
    outside = by >= 0
    new_start = np.where(
        at_start,
        np.where(outside, start - by, start),
        np.where(outside, end + 1, end + 1 + by),
    )
    return new_start, a


def promoters(
    start: np.ndarray,
    width: np.ndarray,
    upstream: np.ndarray,
    downstream: np.ndarray,
    neg: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Promoter windows; `neg` marks ranges whose upstream side is the right end."""
    if (upstream < 0).any() or (downstream < 0).any():
        raise ValueError("'upstream' and 'downstream' must be integers >= 0")
    end = start + width - 1
    pos_start = start - upstream
    if neg is None:
        return pos_start, upstream + downstream
    return np.where(neg, end - downstream + 1, pos_start), upstream + downstream
