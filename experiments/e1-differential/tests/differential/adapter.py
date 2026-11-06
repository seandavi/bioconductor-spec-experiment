"""R call vocabulary -> Python API, and Python value -> canonical JSON.

Add a new R function: usually nothing — an R generic `fooBar(x, ...)` maps to
`x.foo_bar(...)` (kwargs `a.b` -> `a_b`, `fn<-` -> `set_fn`) and properties are
read as attributes. Only add an entry to ADAPTERS when the Python shape differs
(index conventions, non-method functions).
"""

from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np
import polars as pl

import ranges_e1 as r1
from codec import decode


class _Missing:
    """Stand-in for classes the package does not export (yet)."""


GRanges = getattr(r1, "GRanges", _Missing)
Seqinfo = getattr(r1, "Seqinfo", _Missing)
Hits = getattr(r1, "Hits", _Missing)
Rle = getattr(r1, "Rle", _Missing)

# ---- R index <-> Python index ----------------------------------------------


def snake(name: str) -> str:
    name = name.replace(".", "_")
    return re.sub(r"(?<=[a-z0-9])([A-Z])", lambda m: "_" + m.group(1).lower(), name).lower()


def r_subscript(i: Any, n: int) -> Any:
    """R `[` subscript (1-based, negative = exclude, 0 dropped) -> Python subscript."""
    if isinstance(i, list):  # character
        return i
    a = np.asarray(i)
    if np.ma.isMaskedArray(i) and np.ma.getmaskarray(i).any():
        raise IndexError("subscript contains NAs")
    a = np.ma.getdata(a) if np.ma.isMaskedArray(a) else a
    if a.dtype == bool:
        return a
    a = a.astype(np.int64)
    if (a < 0).any():
        if (a > 0).any():
            raise IndexError("only 0's may be mixed with negative subscripts")
        keep = np.ones(n, dtype=bool)
        keep[(-a[a < 0] - 1)[(-a[a < 0] - 1) < n]] = False
        return keep
    a = a[a != 0]
    return a - 1


def to_one_based(v: Any) -> Any:
    """Index-valued Python result -> R's 1-based convention (masked stays NA)."""
    if isinstance(v, np.ma.MaskedArray):
        return np.ma.masked_array(np.ma.getdata(v) + 1, mask=np.ma.getmask(v))
    if isinstance(v, np.ndarray):
        return v + 1
    return v


@dataclass(frozen=True)
class Adapter:
    call: Callable[[list[Any], dict[str, Any]], Any] | None = None
    post: Callable[[Any, list[Any], dict[str, Any]], Any] | None = None


def _subset(args: list[Any], kw: dict[str, Any]) -> Any:
    x, i = args
    return x[r_subscript(i, len(x))]


def _concat(args: list[Any], kw: dict[str, Any]) -> Any:
    return type(args[0]).concat(args)


def _post_select(res: Any, args: list[Any], kw: dict[str, Any]) -> Any:
    return res if kw.get("select", "all") == "all" else to_one_based(res)


def _post_idx(res: Any, args: list[Any], kw: dict[str, Any]) -> Any:
    return to_one_based(res)


ADAPTERS: dict[str, Adapter] = {
    "[": Adapter(call=_subset),
    "c": Adapter(call=_concat),
    "identity": Adapter(call=lambda a, k: a[0]),
    "length": Adapter(call=lambda a, k: len(a[0])),
    "findOverlaps": Adapter(post=_post_select),
    "nearest": Adapter(post=_post_idx),
    "precede": Adapter(post=_post_idx),
    "follow": Adapter(post=_post_idx),
    "order": Adapter(post=_post_idx),
}


# R returns these with the query's names; Python returns plain arrays (ADR-0004).
R_UNNAMED = {"countOverlaps"}


def normalize_r(fn: str | None, value: Any) -> Any:
    """Strip R-only vector names where an ADR says Python does not carry them."""
    if fn in R_UNNAMED and isinstance(value, dict) and value.get("class") == "vector":
        return {**value, "names": None}
    return value


def call_fn(fn: str, args: list[Any], kwargs: dict[str, Any]) -> Any:
    ad = ADAPTERS.get(fn, Adapter())
    kw = {snake(k): v for k, v in kwargs.items()}
    if ad.call is not None:
        res = ad.call(args, kw)
    else:
        name = snake(fn)
        if fn.endswith("<-"):
            name = "set_" + snake(fn[:-2])
            args = [*args, kw.pop("value")] if "value" in kw else args
        target = getattr(args[0], name)
        res = target(*args[1:], **kw) if callable(target) else target
    return ad.post(res, args, kw) if ad.post else res


# ---- evaluation of a case's call tree ---------------------------------------


def _resolve(v: Any, env: dict[str, Any]) -> Any:
    if isinstance(v, str) and v.startswith("$"):
        if v.startswith("$$"):
            return v[1:]
        return env[v[1:]]
    if isinstance(v, dict) and list(v) == ["call"]:
        return run_call(v["call"], env)
    return decode(v)


def run_call(cl: dict[str, Any], env: dict[str, Any]) -> Any:
    args = [_resolve(a, env) for a in cl.get("args", [])]
    kwargs = {k: _resolve(v, env) for k, v in (cl.get("kwargs") or {}).items()}
    return call_fn(cl["fn"], args, kwargs)


def run_python(case: dict[str, Any]) -> Any:
    """Evaluate a case in Python; returns the raw Python value."""
    env = {k: decode(v) for k, v in case.get("inputs", {}).items()}
    if "raw" in case:
        ns = {"np": np, "pl": pl, "r1": r1, **vars(r1), **env}
        return eval(case["raw"]["py"], ns)  # noqa: S307 - trusted case files
    return run_call(case["call"], env)


# ---- Python value -> canonical JSON -----------------------------------------


def _vec_type(a: np.ndarray) -> str:
    return {"b": "lgl", "i": "int", "u": "int", "f": "num"}.get(a.dtype.kind) or "chr"


def enc_vec(values: Any, names: list[str] | None = None) -> dict[str, Any]:
    if isinstance(values, np.ma.MaskedArray):
        typ = _vec_type(values)
        data = [None if m else v for v, m in zip(np.ma.getdata(values).tolist(), np.ma.getmaskarray(values))]
    elif isinstance(values, np.ndarray):
        typ, data = _vec_type(values), values.tolist()
    else:
        data = list(values)
        non_na = [d for d in data if d is not None]
        if non_na and all(isinstance(d, str) for d in non_na) or not non_na:
            typ = "chr"
        elif all(isinstance(d, bool) for d in non_na):
            typ = "lgl"
        elif any(isinstance(d, float) for d in non_na):
            typ = "num"
        else:
            typ = "int"
    return {"class": "vector", "type": typ, "values": data, "names": names}


def _enc_series(s: pl.Series, name: str) -> dict[str, Any]:
    if isinstance(s.dtype, pl.List):
        shift = 1 if name == "revmap" else 0  # revmap holds 0-based member indices
        return {"class": "List", "elements": [[v + shift for v in e] for e in s.to_list()]}
    typ = (
        "lgl" if s.dtype == pl.Boolean
        else "int" if s.dtype.is_integer()
        else "num" if s.dtype.is_float()
        else "chr"
    )
    return {"class": "vector", "type": typ, "values": s.to_list(), "names": None}


def enc_mcols(df: pl.DataFrame | None) -> dict[str, Any] | None:
    if df is None or df.width == 0:
        return None
    return {c: _enc_series(df[c], c) for c in df.columns}


def encode(v: Any) -> Any:
    """Python value -> canonical JSON-able structure (R's 1-based index convention)."""
    if v is None:
        return None
    if isinstance(v, GRanges):
        return v._to_canonical()  # noqa: SLF001
    if isinstance(v, r1.IRanges):
        return {
            "class": "IRanges", "start": v.start.tolist(), "width": v.width.tolist(),
            "names": None if v.names is None else list(v.names), "mcols": enc_mcols(v.mcols),
        }
    if isinstance(v, Seqinfo):
        return v._to_canonical()  # noqa: SLF001
    if isinstance(v, Hits):
        return {
            "class": "Hits", "from": (v.query_hits + 1).tolist(), "to": (v.subject_hits + 1).tolist(),
            "nLnode": v.query_length, "nRnode": v.subject_length,
        }
    if isinstance(v, Rle):
        return {"class": "Rle", "values": v.values.tolist(), "lengths": v.lengths.tolist()}
    if isinstance(v, pl.DataFrame):
        return {"class": "DataFrame", "columns": enc_mcols(v), "nrow": v.height}
    if isinstance(v, dict) and all(isinstance(x, Rle) for x in v.values()) and v:
        return {"class": "RleList", "names": list(v), "elements": [encode(x) for x in v.values()]}
    if isinstance(v, dict):
        return enc_vec(list(v.values()), list(v))
    if isinstance(v, list) and v and all(isinstance(x, (r1.IRanges, GRanges)) for x in v):
        return {"class": "List", "elements": [encode(x) for x in v]}
    if isinstance(v, (bool, int, float, str, np.generic)):
        return enc_vec([v.item() if isinstance(v, np.generic) else v])
    if isinstance(v, (np.ndarray, list, tuple)):
        return enc_vec(v)
    raise TypeError(f"cannot encode {type(v).__name__}")
