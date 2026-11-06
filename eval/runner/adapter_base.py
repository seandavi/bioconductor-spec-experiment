"""Adapter base: binding, call-tree evaluation, projection helpers. stdlib + numpy only.

Contract: ADAPTERS.md. An adapter subclasses AdapterBase and fills in build_* / op_* / project.
This module never imports a port.
"""

from __future__ import annotations

import inspect
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from binding import BindError, bind

OPS_JSON = Path(__file__).resolve().parent.parent / "spec" / "ops.json"
ALL_OPS: dict[str, dict] = {o["op"]: o for o in json.loads(OPS_JSON.read_text())["ops"]}


# ---------------------------------------------------------------- signals

class Unsupported(Exception):
    """The port's public API lacks this op / argument form / input kind."""


class UnknownArgument(Exception):
    """An R argument the adapter cannot map to the port."""

    def __init__(self, name: str):
        super().__init__(name)
        self.name = name


class _Missing:
    def __repr__(self) -> str:
        return "MISSING"

    def __bool__(self) -> bool:
        return False


MISSING = _Missing()


# ---------------------------------------------------------------- atomic argument values

_SPECIAL = {"NaN": math.nan, "Inf": math.inf, "-Inf": -math.inf}


class RVector:
    """An R atomic vector as handed to an adapter. `values` holds JSON-level values (None = NA)."""

    __slots__ = ("type", "values", "names")

    def __init__(self, type: str, values: list, names: list | None = None):
        self.type = type
        self.values = values
        self.names = names

    def __len__(self) -> int:
        return len(self.values)

    def __repr__(self) -> str:
        return f"RVector({self.type}, {self.values!r}, names={self.names!r})"

    def list(self) -> list:
        """Python natives; None = NA; dbl "NaN"/"Inf"/"-Inf" become floats."""
        if self.type == "dbl":
            return [_SPECIAL.get(v, v) if isinstance(v, str) else v for v in self.values]
        return list(self.values)

    def scalar(self):
        """The single element (None = NA). Raises unless length 1."""
        if len(self.values) != 1:
            raise ValueError(f"expected a length-1 R vector, got {len(self.values)}")
        return self.list()[0]

    def np(self, na_as_nan: bool = False) -> np.ndarray:
        """ndarray: int -> int64, dbl -> float64, lgl -> bool, chr -> object.

        Raises ValueError if the vector holds NA; use np_masked() then. A dbl vector may pass
        na_as_nan=True to hand the port NaN for NA (a port with no NA representation should
        instead report the input unsupported).
        """
        vals = self.list()
        if self.type == "dbl" and na_as_nan:
            return np.array([math.nan if v is None else v for v in vals], dtype=np.float64)
        if any(v is None for v in vals):
            raise ValueError(f"NA in {self.type} vector; use np_masked()")
        if self.type == "dbl":
            return np.array(vals, dtype=np.float64)
        dt = {"int": np.int64, "lgl": np.bool_}.get(self.type, object)
        return np.array(vals, dtype=dt)

    def np_masked(self) -> np.ma.MaskedArray:
        """Masked array: NA masked; dbl NaN stays unmasked."""
        vals = self.list()
        mask = [v is None for v in vals]
        fill = {"int": 0, "dbl": 0.0, "lgl": False, "chr": ""}[self.type]
        dt = {"int": np.int64, "dbl": np.float64, "lgl": np.bool_}.get(self.type, object)
        return np.ma.masked_array(np.array([fill if m else v for v, m in zip(vals, mask)], dtype=dt), mask=mask)


def _rvector_from_literal(lit: dict) -> RVector:
    kind = next(k for k in ("int", "dbl", "chr", "lgl") if k in lit)
    return RVector(kind, list(lit[kind]), list(lit["names"]) if lit.get("names") is not None else None)


def decode_arg_literal(a: Any) -> Any:
    """Case-schema atomic argument -> RVector / None (refs and calls are handled by the caller)."""
    if a is None:
        return None
    if isinstance(a, bool):
        return RVector("lgl", [a])
    if isinstance(a, int):
        return RVector("int", [a])
    if isinstance(a, float):
        return RVector("dbl", [a])
    if isinstance(a, str):
        return RVector("chr", [a])
    return _rvector_from_literal(a)


# ---------------------------------------------------------------- projection helpers (PROJECTION.md)

def _py(v: Any, type_: str) -> Any:
    """One scalar -> JSON scalar of an R type; NA (None / masked / nan-in-int) -> None."""
    if v is None or v is np.ma.masked:
        return None
    if type_ == "int":
        if isinstance(v, float):
            if math.isnan(v):
                return None
            if not v.is_integer():
                raise ValueError(f"non-integral value {v!r} in int vector")
        if isinstance(v, (str, bytes)):
            raise ValueError(f"string {v!r} in int vector")
        return int(v)
    if type_ == "dbl":
        f = float(v)
        if math.isnan(f):
            return "NaN"
        if math.isinf(f):
            return "Inf" if f > 0 else "-Inf"
        return f
    if type_ == "lgl":
        if isinstance(v, float) and math.isnan(v):
            return None
        return bool(v)
    if type_ == "chr":
        if isinstance(v, float) and math.isnan(v):
            return None
        return str(v)
    raise ValueError(f"unknown vector type {type_!r}")


def _values_list(arr: Any) -> list:
    if isinstance(arr, np.ma.MaskedArray):
        return [None if m else x for x, m in zip(arr.data.tolist(), np.ma.getmaskarray(arr).tolist())]
    if isinstance(arr, np.ndarray):
        return arr.tolist()
    return list(arr)


def _names(names: Any, n: int) -> list | None:
    if names is None or n == 0:
        return None
    out = [None if (x is None or (isinstance(x, float) and math.isnan(x))) else str(x) for x in _values_list(names)]
    if len(out) != n:
        raise ValueError(f"names length {len(out)} != {n}")
    return out


def infer_type(arr: Any) -> str:
    if isinstance(arr, np.ndarray):
        k = arr.dtype.kind
        if k in "iu":
            return "int"
        if k == "f":
            return "dbl"
        if k == "b":
            return "lgl"
        if k in "USO":
            if k == "O":
                return infer_type(arr.tolist())
            return "chr"
        raise ValueError(f"cannot infer R type from dtype {arr.dtype}")
    vals = [v for v in _values_list(arr) if v is not None and v is not np.ma.masked]
    if not vals:
        raise ValueError("cannot infer R type of an empty / all-NA sequence; use p_vector(type, ...)")
    if all(isinstance(v, (bool, np.bool_)) for v in vals):
        return "lgl"
    if all(isinstance(v, str) for v in vals):
        return "chr"
    if all(isinstance(v, (int, np.integer)) and not isinstance(v, (bool, np.bool_)) for v in vals):
        return "int"
    if all(isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_)) for v in vals):
        return "dbl"
    raise ValueError(f"cannot infer R type from mixed values {vals[:3]!r}")


def p_vector(type: str, values: Any, names: Any = None) -> dict:
    vals = [_py(v, type) for v in _values_list(values)]
    return {"kind": "vector", "type": type, "values": vals, "names": _names(names, len(vals))}


def p_vector_from(arr: Any, names: Any = None) -> dict:
    if isinstance(arr, np.ma.MaskedArray):
        t = infer_type(arr.data)
    else:
        t = infer_type(arr)
    return p_vector(t, arr, names)


def p_null() -> dict:
    return {"kind": "NULL"}


def _ints(v: Any) -> list:
    return [_py(x, "int") for x in _values_list(v)]


def p_mcols(mcols: Any) -> list:
    return [[str(k), col] for k, col in (mcols or [])]


def p_iranges(start: Any, width: Any, names: Any = None, mcols: Any = None) -> dict:
    s, w = _ints(start), _ints(width)
    if len(s) != len(w):
        raise ValueError("start and width differ in length")
    return {"kind": "IRanges", "start": s, "width": w, "names": _names(names, len(s)), "mcols": p_mcols(mcols)}


def p_seqinfo(seqnames: Any, seqlengths: Any = None, is_circular: Any = None, genome: Any = None) -> dict:
    sn = [str(x) for x in _values_list(seqnames)]
    n = len(sn)

    def col(v, t):
        if v is None:
            return [None] * n
        out = [_py(x, t) for x in _values_list(v)]
        if len(out) != n:
            raise ValueError(f"seqinfo column length {len(out)} != {n}")
        return out

    return {"kind": "Seqinfo", "seqnames": sn, "seqlengths": col(seqlengths, "int"),
            "isCircular": col(is_circular, "lgl"), "genome": col(genome, "chr")}


def p_granges(seqnames: Any, start: Any, width: Any, strand: Any, names: Any = None, mcols: Any = None,
              seqinfo: dict | None = None) -> dict:
    s, w = _ints(start), _ints(width)
    sn = [str(x) for x in _values_list(seqnames)]
    st = [str(x) for x in _values_list(strand)]
    if not (len(s) == len(w) == len(sn) == len(st)):
        raise ValueError("seqnames/start/width/strand differ in length")
    if seqinfo is None or seqinfo.get("kind") != "Seqinfo":
        raise ValueError("p_granges needs a p_seqinfo(...) projection")
    return {"kind": "GRanges", "seqnames": sn, "start": s, "width": w, "strand": st,
            "names": _names(names, len(s)), "mcols": p_mcols(mcols), "seqinfo": seqinfo}


def p_factor(levels: Any, codes: Any, base: int = 1) -> dict:
    """codes are 0-based when base=0 (negative = NA), else 1-based; the projection is 1-based."""
    out = []
    for c in _values_list(codes):
        if c is None or c is np.ma.masked or (isinstance(c, float) and math.isnan(c)):
            out.append(None)
        elif base == 0 and int(c) < 0:
            out.append(None)
        else:
            out.append(int(c) - base + 1)
    return {"kind": "Factor", "levels": [str(x) for x in _values_list(levels)], "codes": out}


def p_integerlist(elements: Any, names: Any = None, shift: int = 0) -> dict:
    """`shift` is added to every non-NA element (use it for index-valued lists such as revmap)."""
    def one(x: Any) -> Any:
        v = _py(x, "int")
        return None if v is None else v + shift

    els = [[one(x) for x in _values_list(e)] for e in elements]
    return {"kind": "IntegerList", "elements": els, "names": _names(names, len(els))}


def p_dataframe(nrow: int, columns: Any, rownames: Any = None) -> dict:
    return {"kind": "DataFrame", "nrow": int(nrow), "rownames": _names(rownames, int(nrow)),
            "columns": p_mcols(columns)}


def p_hits(from_: Any, to: Any, n_left: int, n_right: int) -> dict:
    """Pairs in the port's emitted order; 1-based indices in."""
    f, t = _ints(from_), _ints(to)
    if len(f) != len(t):
        raise ValueError("from and to differ in length")
    return {"kind": "Hits", "from": f, "to": t, "nLnode": int(n_left), "nRnode": int(n_right)}


def _run_key(v: Any):
    if v is None:
        return ("NA",)
    if isinstance(v, float) and math.isnan(v):
        return ("NaN",)
    return ("v", v)


def _rle_body(values: Any, lengths: Any = None, type: str | None = None) -> dict:
    if isinstance(values, np.ma.MaskedArray):
        t = type or infer_type(values.data)
    elif isinstance(values, np.ndarray) and values.dtype.kind == "O":
        t = type or infer_type(values)
    else:
        t = type or infer_type(values)
    if t not in ("int", "dbl", "lgl"):
        raise ValueError(f"Rle type must be int/dbl/lgl, got {t}")
    vals = _values_list(values)
    lens = [1] * len(vals) if lengths is None else [int(x) for x in _values_list(lengths)]
    if len(lens) != len(vals):
        raise ValueError("Rle values and lengths differ in length")
    py = [_py(v, t) for v in vals]
    out_v: list = []
    out_l: list = []
    for v, n in zip(py, lens):
        if n <= 0:
            continue
        if out_v and _run_key(out_v[-1]) == _run_key(v):
            out_l[-1] += n
        else:
            out_v.append(v)
            out_l.append(n)
    return {"type": t, "values": out_v, "lengths": out_l}


def p_rle(values: Any, lengths: Any = None, type: str | None = None) -> dict:
    """Dense (values only) or run form; normalized to maximal runs."""
    return {"kind": "Rle", **_rle_body(values, lengths, type)}


def p_rlelist(names: Any, rles: Any) -> dict:
    """rles: Rle projections, or (values, lengths) tuples, or dense arrays."""
    els = []
    for r in rles:
        if isinstance(r, dict):
            els.append(_rle_body(r["values"], r["lengths"], r["type"]))
        elif isinstance(r, tuple):
            els.append(_rle_body(*r))
        else:
            els.append(_rle_body(r))
    nm = None if names is None else [None if x is None else str(x) for x in _values_list(names)]
    if nm is not None and len(nm) == 0:
        nm = None
    return {"kind": "RleList", "names": nm, "elements": els}


# ---------------------------------------------------------------- the base class

class AdapterBase:
    """Subclass per port. See ADAPTERS.md."""

    port: str = ""
    port_packages: tuple[str, ...] = ()
    unsupported_ops: dict[str, str] = {}

    # ---- inputs ------------------------------------------------------------
    def build_vector(self, p: dict) -> Any:
        return RVector(p["type"], list(p["values"]), list(p["names"]) if p["names"] is not None else None)

    def build_input(self, p: dict) -> Any:
        kind = p["kind"]
        if kind == "NULL":
            return None
        if kind == "vector":
            return self.build_vector(p)
        fn = getattr(self, f"build_{kind}", None)
        if fn is None:
            raise Unsupported(f"input:{kind}")
        return fn(p)

    def build_inputs(self, case: dict) -> dict:
        return {name: self.build_input(p) for name, p in case["inputs"].items()}

    # ---- calls ---------------------------------------------------------------
    def eval_arg(self, a: Any, env: dict) -> Any:
        if isinstance(a, dict):
            if "$ref" in a and len(a) == 1:
                return env[a["$ref"]]
            if "$call" in a and len(a) == 1:
                return self.call(a["$call"], env)
        return decode_arg_literal(a)

    def call(self, call: dict, env: dict) -> Any:
        """Evaluate one Call node: args innermost first, bind to R names, dispatch to op_<method>."""
        op = call["op"]
        entry = ALL_OPS[op]
        args = [self.eval_arg(a, env) for a in call["args"]]
        kwargs = {k: self.eval_arg(v, env) for k, v in call["kwargs"].items()}
        if op in self.unsupported_ops:
            raise Unsupported(f"op:{op}: {self.unsupported_ops[op]}")
        bound, extra_pos, extra_named = bind(entry, args, kwargs)  # BindError: R would reject it too
        kw = {k.replace(".", "_"): v for k, v in bound.items()}
        if extra_named:
            kw["_mcols"] = extra_named
        fn = getattr(self, f"op_{entry['method']}", None)
        if fn is None:
            raise AttributeError(f"adapter {type(self).__name__} defines neither op_{entry['method']} nor unsupported_ops[{op!r}]")
        sig = inspect.signature(fn)
        accepts_kw = any(p.kind is p.VAR_KEYWORD for p in sig.parameters.values())
        if not accepts_kw:
            for k in kw:
                if k not in sig.parameters:
                    raise UnknownArgument(k)
        if extra_pos and not any(p.kind is p.VAR_POSITIONAL for p in sig.parameters.values()):
            raise Unsupported(f"op:{op}: extra positional objects")
        return fn(*extra_pos, **kw)

    def evaluate(self, case: dict, inputs: dict) -> Any:
        return self.call(case["call"], inputs)

    # ---- results --------------------------------------------------------------
    def project(self, value: Any) -> dict:
        raise NotImplementedError("adapter must implement project(value)")

    def _project(self, value: Any) -> dict:
        """RVector and None are projected here; everything else by the adapter."""
        if value is None:
            return p_null()
        if isinstance(value, RVector):
            return p_vector(value.type, value.values, value.names)
        return self.project(value)

    # ---- one case (the runner classifies exceptions raised from these stages) ------
    def run_case(self, case: dict, stage) -> dict:
        """Default pipeline. `stage(phase, fn, *a)` runs fn under phase-tagged classification.

        Returns {"result": projection}. Self-test adapters override this.
        """
        inputs = stage("build", self.build_inputs, case)
        value = stage("call", self.evaluate, case, inputs)
        result = stage("project", self._project, value)
        return {"result": result}
