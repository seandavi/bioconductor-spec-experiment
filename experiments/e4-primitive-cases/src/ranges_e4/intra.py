"""Intra-range kernels over (start, width) int64 arrays.

``reverse`` is a per-range bool mask: True means the range is on the minus
strand, which swaps the roles of start and end (GRanges passes its strand
mask; IRanges passes all-False).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ._arrays import as_int


def arg_vector(value: Any, n: int, name: str) -> np.ndarray:
    """Integer argument: length 1 recycles; otherwise must have length ``n``."""
    arr = as_int(value, name)
    if len(arr) == 1:
        return np.repeat(arr, n)
    if len(arr) != n:
        raise ValueError(f"'{name}' has length {len(arr)}, expected 1 or {n}")
    return arr


def flag_vector(value: Any, n: int, name: str) -> np.ndarray:
    arr = np.atleast_1d(np.asarray(value))
    if arr.dtype != bool:
        raise TypeError(f"'{name}' must be logical")
    if len(arr) == 1:
        return np.repeat(arr, n)
    if len(arr) != n:
        raise ValueError(f"'{name}' has length {len(arr)}, expected 1 or {n}")
    return arr


def str_vector(value: Any, n: int, name: str) -> np.ndarray:
    arr = np.atleast_1d(np.asarray(value, dtype=object))
    if len(arr) == 1:
        return np.repeat(arr, n)
    if len(arr) != n:
        raise ValueError(f"'{name}' has length {len(arr)}, expected 1 or {n}")
    return arr


def shift(start: np.ndarray, width: np.ndarray, by: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return start + by, width


def _resolve(v: np.ndarray, refwidth: np.ndarray) -> np.ndarray:
    """Negative offsets count back from the range end: -1 is the last position."""
    return np.where(v < 0, refwidth + v + 1, v)


def narrow(
    start: np.ndarray,
    width: np.ndarray,
    s: np.ndarray | None,
    e: np.ndarray | None,
    w: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    """R ``narrow`` (``solveUserSEW`` with ``allow.nonnarrowing=FALSE``)."""
    if s is None and e is None:
        raise ValueError("'start' or 'end' must be supplied")
    rw = width
    if s is not None and e is not None:
        if w is not None:
            raise ValueError("'width' must not be supplied with both 'start' and 'end'")
        ns, ne = _resolve(s, rw), _resolve(e, rw)
    elif s is not None:
        ns = _resolve(s, rw)
        ne = rw if w is None else ns + w - 1
    else:
        assert e is not None
        ne = _resolve(e, rw)
        ns = np.ones_like(rw) if w is None else ne - w + 1
    if ((ns < 1) | (ns > rw + 1) | (ne > rw) | (ne < 0) | (ne < ns - 1)).any():
        raise ValueError("narrowing would extend beyond the original range")
    return start + ns - 1, ne - ns + 1


def resize(
    start: np.ndarray, width: np.ndarray, new_width: np.ndarray, fix: np.ndarray, reverse: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    if (new_width < 0).any():
        raise ValueError("'width' values must be non-negative")
    codes = np.full(len(fix), -1)
    for name, code in (("start", 0), ("end", 1), ("center", 2)):
        codes[fix == name] = code
    if (codes < 0).any():
        raise ValueError("'fix' must be 'start', 'end' or 'center'")
    codes = np.where(reverse & (codes < 2), 1 - codes, codes)
    end = start + width - 1
    new_start = np.where(
        codes == 0, start, np.where(codes == 1, end - new_width + 1, start + (width - new_width) // 2)
    )
    return new_start, new_width


def flank(
    start: np.ndarray, width: np.ndarray, w: np.ndarray, at_start: np.ndarray, both: bool, reverse: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    end = start + width - 1
    at_start = at_start ^ reverse
    if both:
        aw = np.abs(w)
        return np.where(at_start, start - aw, end - aw + 1), 2 * aw
    ns = np.where(
        at_start,
        np.where(w >= 0, start - w, start),
        np.where(w >= 0, end + 1, end + w + 1),
    )
    return ns, np.abs(w)


def promoters(
    start: np.ndarray, width: np.ndarray, up: np.ndarray, down: np.ndarray, reverse: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    if (up < 0).any() or (down < 0).any():
        raise ValueError("'upstream' and 'downstream' must be non-negative")
    end = start + width - 1
    return np.where(reverse, end - down + 1, start - up), up + down
