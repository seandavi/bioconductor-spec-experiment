"""Argument validation mirroring R's `isTRUEorFALSE`, `isSingleNumber` and `match.arg` (internal)."""

from __future__ import annotations

import functools
import inspect
from collections.abc import Callable, Sequence
from typing import Any, TypeVar

import numpy as np

from ._meta import INT_MAX

F = TypeVar("F", bound=Callable[..., Any])


def _single(v: Any) -> Any:
    """Unwrap a length-1 list/array to its element; anything else of non-scalar shape is returned as is."""
    if isinstance(v, np.ma.MaskedArray):
        if v.size != 1:
            return v
        return None if np.ma.getmaskarray(v).ravel()[0] else np.ma.getdata(v).ravel()[0]
    if isinstance(v, (list, tuple, np.ndarray)) and np.size(v) == 1:
        return np.asarray(v).ravel()[0]
    return v


def flag(v: Any, what: str) -> bool:
    """TRUE / FALSE only (length-1 arrays allowed); NA and non-logicals are errors."""
    x = _single(v)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    raise ValueError(f"'{what}' must be TRUE or FALSE")


def logical_vector(v: Any, what: str) -> np.ndarray:
    """A logical vector without NAs (masked entries included) as a bool ndarray."""
    if isinstance(v, np.ma.MaskedArray) and np.ma.getmaskarray(v).any():
        raise ValueError(f"'{what}' must be logical without NA's")
    a = np.atleast_1d(np.ma.getdata(v) if isinstance(v, np.ma.MaskedArray) else np.asarray(v))
    if a.dtype.kind != "b":
        raise ValueError(f"'{what}' must be logical without NA's")
    return a


def scalar_int(v: Any, what: str, *, allow_na: bool = False) -> int | None:
    """A single number as int (doubles truncate like as.integer); NA only when allowed."""
    x = _single(v)
    if x is None:
        if allow_na:
            return None
        raise ValueError(f"'{what}' must be a single integer")
    if isinstance(x, (bool, np.bool_, str, bytes)) or not isinstance(x, (int, float, np.integer, np.floating)):
        raise ValueError(f"'{what}' must be a single integer")
    if isinstance(x, (float, np.floating)):
        if x != x:
            if allow_na:
                return None
            raise ValueError(f"'{what}' must be a single integer")
        x = int(np.trunc(x))
    x = int(x)
    if abs(x) > INT_MAX:
        raise ValueError(f"'{what}' is outside the integer range")
    return x


def match_arg(v: Any, choices: Sequence[str], what: str = "arg") -> str:
    """R's match.arg: exact match, else a unique prefix; a single string only."""
    x = _single(v)
    if not isinstance(x, (str, np.str_)):
        raise ValueError(f"'{what}' must be a character string")
    if x in choices:
        return str(x)
    hits = [c for c in choices if c.startswith(x)] if x else []
    if len(hits) == 1:
        return hits[0]
    raise ValueError(f"'{what}' should be one of {', '.join(choices)}")


def flags(*names: str) -> Callable[[F], F]:
    """Decorator: validate the named (keyword or positional) arguments as TRUE/FALSE flags."""

    def deco(fn: F) -> F:
        sig = inspect.signature(fn)

        @functools.wraps(fn)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            bound = sig.bind(*args, **kwargs)  # unexpected keywords raise TypeError, as R's "unused argument"
            for n in names:
                if n in bound.arguments:
                    bound.arguments[n] = flag(bound.arguments[n], n.replace("_", "."))
            return fn(*bound.args, **bound.kwargs)

        return wrapper  # type: ignore[return-value]

    return deco
