"""Shared helpers: integer coercion, start/end/width solving, indexing, mcols."""

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import polars as pl

type IntArray = npt.NDArray[np.int64]
type BoolArray = npt.NDArray[np.bool_]

NO_HIT = -1
"""Missing index (R `NA`) in positional results such as `select="first"`."""


def empty_int() -> IntArray:
    return np.empty(0, dtype=np.int64)


def as_int_array(x: Any, name: str = "value") -> IntArray:
    """Coerce a scalar or sequence of integers to a 1-d int64 array."""
    a = np.asarray(x)
    if a.ndim == 0:
        a = a.reshape(1)
    if a.ndim != 1:
        raise ValueError(f"{name} must be 1-dimensional")
    kind = a.dtype.kind
    if kind in "iu":
        return a.astype(np.int64, copy=False)
    if kind == "f":
        if not np.all(np.isfinite(a)) or np.any(a != np.floor(a)):
            raise ValueError(f"{name} must contain whole numbers")
        return a.astype(np.int64)
    raise ValueError(f"{name} must be integer, got dtype {a.dtype}")


def per_range(x: Any, n: int, name: str) -> IntArray:
    """Scalar or length-n integer argument broadcast to length n."""
    a = as_int_array(x, name)
    if len(a) == n:
        return a
    if len(a) == 1:
        return np.repeat(a, n)
    raise ValueError(f"{name} must have length 1 or {n}, got {len(a)}")


def common_length(arrays: Sequence[IntArray], names: Sequence[str]) -> int:
    lens = {len(a) for a in arrays}
    others = lens - {1}
    if len(others) > 1:
        raise ValueError(
            "incompatible lengths: "
            + ", ".join(f"{nm}={len(a)}" for nm, a in zip(names, arrays))
        )
    return others.pop() if others else 1


def solve_sew(
    start: Any = None, end: Any = None, width: Any = None
) -> tuple[IntArray, IntArray]:
    """Resolve any two of start/end/width to (start, width); width must be >= 0."""
    given = [(nm, v) for nm, v in (("start", start), ("end", end), ("width", width)) if v is not None]
    if not given:
        return empty_int(), empty_int()
    if len(given) == 1:
        raise ValueError("at least two of start, end, width are required")
    arrs = {nm: as_int_array(v, nm) for nm, v in given}
    n = common_length(list(arrs.values()), list(arrs))
    arrs = {nm: (np.repeat(a, n) if len(a) != n else a) for nm, a in arrs.items()}
    if "start" in arrs and "end" in arrs:
        s = arrs["start"]
        w = arrs["end"] - s + 1
        if "width" in arrs and np.any(arrs["width"] != w):
            raise ValueError("start, end and width are inconsistent")
    elif "start" in arrs:
        s, w = arrs["start"], arrs["width"]
    else:
        w = arrs["width"]
        s = arrs["end"] - w + 1
    if np.any(w < 0):
        raise ValueError("ranges must have non-negative width (end >= start - 1)")
    return s, w


def resolve_index(i: Any, n: int, names: Sequence[str] | None) -> IntArray:
    """Turn int / slice / int array / bool mask / name(s) into positions (0-based)."""
    if isinstance(i, (int, np.integer)) and not isinstance(i, (bool, np.bool_)):
        j = int(i)
        if j < -n or j >= n:
            raise IndexError(f"index {j} out of range for length {n}")
        return np.array([j % n], dtype=np.int64)
    if isinstance(i, slice):
        return np.arange(n, dtype=np.int64)[i]
    if isinstance(i, str):
        i = [i]
    if isinstance(i, (list, tuple)) and i and all(isinstance(v, str) for v in i):
        if names is None:
            raise KeyError("object has no names")
        first: dict[str, int] = {}
        for pos in range(len(names) - 1, -1, -1):
            first[names[pos]] = pos
        try:
            return np.array([first[v] for v in i], dtype=np.int64)
        except KeyError as err:
            raise KeyError(f"name not found: {err.args[0]!r}") from None
    a = np.asarray(i)
    if a.dtype.kind == "b":
        if a.shape != (n,):
            raise ValueError(f"boolean mask must have length {n}, got shape {a.shape}")
        return np.flatnonzero(a).astype(np.int64)
    if a.size == 0:
        return empty_int()
    if a.dtype.kind in "iu":
        a = a.astype(np.int64, copy=False).ravel()
        if np.any(a < -n) or np.any(a >= n):
            raise IndexError(f"index out of range for length {n}")
        return np.where(a < 0, a + n, a)
    raise TypeError(f"unsupported index type: {type(i).__name__}")


# --- mcols -----------------------------------------------------------------


def normalize_mcols(mcols: pl.DataFrame | Mapping[str, Any] | None, n: int) -> pl.DataFrame:
    if mcols is None:
        return pl.DataFrame()
    df = mcols.clone() if isinstance(mcols, pl.DataFrame) else pl.DataFrame(mcols)
    if df.width and df.height != n:
        raise ValueError(f"mcols has {df.height} rows, expected {n}")
    return df


def take_mcols(mcols: pl.DataFrame, idx: IntArray) -> pl.DataFrame:
    if mcols.width == 0:
        return mcols
    return mcols[idx]


def concat_mcols(frames: Sequence[pl.DataFrame]) -> pl.DataFrame:
    """Vertical concat; every input must have identical column names (R rbind rule)."""
    cols = [tuple(f.columns) for f in frames]
    if any(c != cols[0] for c in cols):
        raise ValueError(f"cannot combine: mcols columns differ: {cols}")
    if not cols[0]:
        return pl.DataFrame()
    return pl.concat(list(frames), how="vertical")


def concat_names(
    names: Sequence[list[str] | None], lengths: Sequence[int]
) -> list[str] | None:
    """Names are kept if any input has them; others contribute empty strings."""
    if all(nm is None for nm in names):
        return None
    out: list[str] = []
    for nm, n in zip(names, lengths):
        out.extend(nm if nm is not None else [""] * n)
    return out
