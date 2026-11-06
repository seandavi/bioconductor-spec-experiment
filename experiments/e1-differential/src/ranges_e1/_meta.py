"""Shared helpers for names, mcols and R-style indexing (internal)."""

from __future__ import annotations

from collections.abc import Sequence
import warnings
from typing import Any

import numpy as np
import polars as pl

Names = tuple[str | None, ...]
INT_MAX = 2**31 - 1


def freeze(a: Any, dtype: Any) -> np.ndarray:
    """A read-only array nobody can re-enable for writing.

    It is backed by an immutable `bytes` buffer; arrays already so backed are shared, anything
    else (including read-only views of someone's writable array) is copied.
    """
    arr = np.asarray(a, dtype=dtype)
    if isinstance(arr.base, bytes) and not arr.flags.writeable and arr.flags.c_contiguous:
        return arr
    return np.frombuffer(np.ascontiguousarray(arr, dtype=dtype).tobytes(), dtype=dtype)


def readonly_view(a: np.ndarray) -> np.ndarray:
    """Accessor result: a read-only view, so the stored array cannot be re-enabled for writing."""
    v = a.view()
    v.flags.writeable = False
    return v


def as_int_array(x: Any, what: str = "value") -> np.ndarray:
    """Coerce a numeric vector to int64 like R's as.integer: doubles truncate; NA, logical or text is an error."""
    if np.ma.isMaskedArray(x):
        if np.ma.getmaskarray(x).any():
            raise ValueError(f"{what} cannot contain NAs")
        x = np.ma.getdata(x)
    a = np.asarray(x)
    if a.dtype == object:
        if any(v is None for v in a.ravel()):
            raise ValueError(f"{what} cannot contain NAs")
        a = a.astype(float) if any(isinstance(v, float) for v in a.ravel()) else a.astype(np.int64)
    a = np.atleast_1d(a)
    if a.ndim != 1:
        raise ValueError(f"{what} must be 1-dimensional")
    if a.dtype.kind == "f":
        if np.isnan(a).any() or (np.abs(a) > INT_MAX).any():
            raise ValueError(f"{what} cannot contain NAs")
        a = np.trunc(a).astype(np.int64)
    elif a.dtype.kind not in "iu":
        raise ValueError(f"{what} must be numeric vectors")
    return a.astype(np.int64, copy=False)


def recycle(arrays: Sequence[np.ndarray | None]) -> list[np.ndarray | None]:
    """Recycle non-None arrays to the longest length, as R does."""
    lens = [len(a) for a in arrays if a is not None]
    n = 0 if 0 in lens else max(lens, default=0)
    out: list[np.ndarray | None] = []
    for a in arrays:
        if a is None:
            out.append(None)
        elif len(a) == n:
            out.append(a)
        elif len(a) == 0:
            raise ValueError("cannot recycle a zero-length argument")
        else:
            out.append(np.resize(a, n))
    return out


def recycle_to(a: np.ndarray, n: int, what: str) -> np.ndarray:
    """Recycle `a` to length `n`; zero-length only when n == 0."""
    if len(a) == n:
        return a
    if n == 0:
        return a[:0]
    if len(a) == 0:
        raise ValueError(f"'{what}' has length 0 but 'x' has length {n}")
    return np.resize(a, n)


def recycle_arg(
    a: np.ndarray, n: int, what: str, *, longer: str = "error", warn: bool = True
) -> np.ndarray:
    """Recycle an argument against a length-`n` object, R style.

    Zero-length args are an error unless n == 0. `longer`: "error" or "truncate".
    A non-multiple length warns (`warn=False` to stay silent).
    """
    m = len(a)
    if m == n:
        return a
    if n == 0:
        if m > 1 and longer == "error":
            raise ValueError(f"'{what}' is longer than the object")
        return a[:0]
    if m == 0:
        raise ValueError(f"'{what}' has no elements")
    if m > n and longer == "error":
        raise ValueError(f"'{what}' is longer than the object")
    if warn and (n % m if m < n else m % n):
        warnings.warn(f"'{what}' length is not a multiple of the object length", UserWarning, stacklevel=3)
    return np.resize(a, n)


def _name_text(v: Any) -> str | None:
    """R's as.character for one element of a names vector (None = NA)."""
    if v is None or v is np.ma.masked:
        return None
    if isinstance(v, (bool, np.bool_)):
        return "TRUE" if v else "FALSE"
    if isinstance(v, (float, np.floating)):
        if v != v:
            return None
        return str(int(v)) if float(v).is_integer() and abs(v) < 1e15 else f"{float(v):.15g}"
    return str(v)


def check_names(names: Any, n: int) -> Names | None:
    if names is None:
        return None
    seq = [_name_text(v) for v in list(names)]
    if len(seq) > n:
        raise ValueError(f"attempt to set too many names ({len(seq)}) on object of length {n}")
    seq.extend([None] * (n - len(seq)))
    return tuple(seq)


def check_mcols(mcols: Any, n: int) -> pl.DataFrame | None:
    if mcols is None:
        return None
    df = pl.DataFrame(mcols) if not isinstance(mcols, pl.DataFrame) else mcols.clone()  # detach from caller's frame
    if df.width == 0:
        return None
    if df.height == n:
        return df
    if df.height == 0 or df.height > n:
        raise ValueError(f"mcols has {df.height} rows but object has length {n}")
    if n % df.height:
        warnings.warn(
            f"You supplied a metadata column of length {df.height} to set on an object of length {n}. "
            "However please note that the latter is not a multiple of the former.",
            UserWarning, stacklevel=4,
        )
    return df[(np.arange(n) % df.height).tolist()]  # R recycles short metadata rows


def take_names(names: Names | None, pos: np.ndarray) -> Names | None:
    if names is None:
        return None
    return tuple(names[i] for i in pos)


def take_mcols(mcols: pl.DataFrame | None, pos: np.ndarray) -> pl.DataFrame | None:
    if mcols is None:
        return None
    return mcols[pos.tolist()] if len(pos) else mcols.clear()


def concat_names(parts: Sequence[tuple[Names | None, int]]) -> Names | None:
    """Concatenate names; unnamed parts contribute "" (R behaviour)."""
    if all(nm is None for nm, _ in parts):
        return None
    out: list[str | None] = []
    for nm, n in parts:
        out.extend(nm if nm is not None else [""] * n)
    return tuple(out)


def concat_mcols(parts: Sequence[tuple[pl.DataFrame | None, int]]) -> pl.DataFrame | None:
    """Row-bind mcols; columns = union in first-seen order, NA-filled."""
    cols: dict[str, pl.DataType] = {}
    for df, _ in parts:
        if df is not None:
            for name, dt in df.schema.items():
                cols.setdefault(name, dt)
    if not cols:
        return None
    series: list[pl.Series] = []
    for name, dt in cols.items():
        pieces = [
            df[name] if df is not None and name in df.columns else pl.Series(name, [None] * n, dtype=dt)
            for df, n in parts
        ]
        if any(p.dtype == pl.String for p in pieces):  # R coerces to character with as.character()
            pieces = [
                p if p.dtype == pl.String else pl.Series(name, [_name_text(v) for v in p.to_list()], dtype=pl.String)
                for p in pieces
            ]
        series.append(pl.concat([p.to_frame() for p in pieces], how="vertical_relaxed")[name])
    return pl.DataFrame(series)


def resolve_index(i: Any, n: int, names: Names | None) -> np.ndarray:
    """Python-style subscript -> int64 positions (0-based, validated)."""
    if isinstance(i, slice):
        return np.arange(n, dtype=np.int64)[i]
    if isinstance(i, (int, np.integer)) and not isinstance(i, (bool, np.bool_)):
        i = [int(i)]
    if isinstance(i, str):
        i = [i]
    if np.ma.isMaskedArray(i):
        raise IndexError("subscript contains NAs")
    a = np.asarray(i)
    if a.size == 0:
        if a.dtype.kind in "US" and names is None:
            raise IndexError("cannot subset by name when names are None")
        return np.empty(0, dtype=np.int64)
    if a.dtype.kind == "b":
        if a.ndim != 1:
            raise IndexError("logical subscript must be 1-dimensional")
        if len(a) > n and a[n:].any():
            raise IndexError("subscript is a logical vector with out-of-bounds TRUE values")
        mask = np.resize(a, n) if 0 < len(a) < n else a[:n]  # R recycles short logical subscripts
        return np.flatnonzero(mask).astype(np.int64)
    if a.dtype.kind in "US" or a.dtype == object and any(isinstance(v, str) for v in a.ravel()):
        if names is None:
            raise IndexError("cannot subset by name when names are None")
        lookup: dict[str, int] = {}
        for k, nm in enumerate(names):
            if nm is not None:
                lookup.setdefault(nm, k)
        try:
            return np.array([lookup[str(v)] for v in a.ravel()], dtype=np.int64)
        except KeyError as e:
            raise IndexError(f"subscript contains invalid names: {e}") from None
    if a.dtype.kind == "f":
        if np.isnan(a).any():
            raise IndexError("subscript contains NAs")
        a = a.astype(np.int64)
    if a.dtype.kind not in "iu":
        raise IndexError("unsupported subscript type")
    a = a.astype(np.int64).ravel()
    if ((a >= n) | (a < -n)).any():
        raise IndexError("subscript contains out-of-bounds indices")
    return np.where(a < 0, a + n, a)


def revmap_frame(members: Sequence[np.ndarray]) -> pl.DataFrame:
    """One-column mcols holding 0-based member positions (R's `revmap`)."""
    col = pl.Series("revmap", [m.tolist() for m in members], dtype=pl.List(pl.Int64))
    return pl.DataFrame([col])
