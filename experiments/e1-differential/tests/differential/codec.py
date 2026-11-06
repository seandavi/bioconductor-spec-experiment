"""Case-literal decoding (JSON/YAML -> Python objects) and canonical-value comparison."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import polars as pl

import ranges_e1 as r1

_TYPED = {"int": np.int64, "num": np.float64, "lgl": np.bool_}


def decode_vec(v: list[Any], typ: str | None = None) -> Any:
    """Array literal -> ndarray (numeric/logical, masked if NA) or list[str|None]."""
    if typ == "chr" or (typ is None and any(isinstance(e, str) for e in v)):
        return [None if e is None else str(e) for e in v]
    if not v:
        return np.empty(0, dtype=_TYPED.get(typ or "int"))
    has_na = any(e is None for e in v)
    non_na = [e for e in v if e is not None]
    if typ is None:
        if all(isinstance(e, bool) for e in non_na) and non_na:
            typ = "lgl"
        elif any(isinstance(e, float) for e in non_na):
            typ = "num"
        else:
            typ = "int"
    dt = _TYPED[typ]
    fill = 0
    arr = np.array([fill if e is None else e for e in v], dtype=dt)
    if has_na:
        return np.ma.masked_array(arr, mask=[e is None for e in v])
    return arr


def _series(name: str, values: Any) -> pl.Series:
    if isinstance(values, np.ma.MaskedArray):
        return pl.Series(name, values.filled(0).tolist()).zip_with(
            pl.Series(~np.ma.getmaskarray(values)), pl.Series([None] * len(values))
        )
    if isinstance(values, np.ndarray):
        return pl.Series(name, values.tolist())
    return pl.Series(name, values)


def _g(d: dict[str, Any], k: str) -> Any:
    return decode(d[k]) if d.get(k) is not None else None


def decode(v: Any) -> Any:
    """Decode one literal (see tests/differential/README.md)."""
    if isinstance(v, dict):
        k = next(iter(v))
        body = v[k]
        if k in _TYPED or k == "chr":
            return decode_vec(body, k)
        if k == "IRanges":
            return r1.IRanges(
                start=_g(body, "start"), end=_g(body, "end"), width=_g(body, "width"),
                names=_g(body, "names"), mcols=_mcols_lit(body.get("mcols")),
            )
        if k == "GRanges":
            return r1.GRanges(
                _g(body, "seqnames"), start=_g(body, "start"), end=_g(body, "end"),
                width=_g(body, "width"), strand=_g(body, "strand"), names=_g(body, "names"),
                mcols=_mcols_lit(body.get("mcols")),
                seqinfo=decode({"Seqinfo": body["seqinfo"]}) if body.get("seqinfo") else None,
            )
        if k == "Seqinfo":
            return r1.Seqinfo(
                seqnames=_g(body, "seqnames"), seqlengths=_g(body, "seqlengths"),
                is_circular=_g(body, "is_circular"), genome=_g(body, "genome"),
            )
        if k == "DataFrame":
            return _mcols_lit(body)
        raise ValueError(f"unknown literal: {k}")
    if isinstance(v, list):
        return decode_vec(v)
    return v


def _mcols_lit(cols: dict[str, Any] | None) -> pl.DataFrame | None:
    if not cols:
        return None
    return pl.DataFrame([_series(name, decode(v)) for name, v in cols.items()])


# ---- comparison -------------------------------------------------------------

_NUM = {"int", "num"}


def _isnum(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def first_diff(r: Any, p: Any, path: str = "$") -> str | None:
    """First difference between two canonical values, or None if equal.

    int/num compare by value (|d| <= 1e-9); everything else strictly.
    Errors compare equal if both sides are errors (messages ignored).
    """
    if isinstance(r, dict) and isinstance(p, dict):
        if "error" in r or "error" in p:
            return None if ("error" in r and "error" in p) else f"{path}: error on one side only"
        if list(r) != list(p) and (set(r) != set(p) or path.endswith("mcols")):
            return f"{path}: keys {list(r)} != {list(p)}"
        for key in r:
            if key == "type" and r[key] in _NUM and p.get(key) in _NUM:
                continue
            d = first_diff(r[key], p[key], f"{path}.{key}")
            if d:
                return d
        return None
    if isinstance(r, list) and isinstance(p, list):
        if len(r) != len(p):
            return f"{path}: length {len(r)} != {len(p)}"
        for i, (a, b) in enumerate(zip(r, p)):
            d = first_diff(a, b, f"{path}[{i}]")
            if d:
                return d
        return None
    if _isnum(r) and _isnum(p):
        if r == p or (math.isfinite(r) and math.isfinite(p) and abs(r - p) <= 1e-9):
            return None
        return f"{path}: {r!r} != {p!r}"
    if r == p and type(r) is type(p):
        return None
    return f"{path}: {r!r} != {p!r}"
