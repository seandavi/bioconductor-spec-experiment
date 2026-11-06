"""Shared converter pieces: unconvertible signalling, literal rules, hoisting, constructor trees."""

from __future__ import annotations

import math
from typing import Any

from ..common import ops_in_call, load_ops

INT32 = 2**31 - 1


class Unconvertible(Exception):
    def __init__(self, reason: str, detail: str = ""):
        super().__init__(reason, detail)
        self.reason = reason
        self.detail = detail


def is_int32(v: Any) -> bool:
    return isinstance(v, int) and not isinstance(v, bool) and -INT32 <= v <= INT32


def dbl_value(v: Any) -> Any:
    """Python number -> element of a dbl vector literal."""
    if v is None:
        return None
    if isinstance(v, float):
        if math.isnan(v):
            return "NaN"
        if math.isinf(v):
            return "Inf" if v > 0 else "-Inf"
        return v
    return float(v)


def vec(kind: str, values: list, names: list | None = None) -> dict:
    d: dict = {kind: values}
    if names is not None:
        d["names"] = names
    return d


def untyped_vector(values: list, *, empty: str | None) -> dict:
    """R `unlist(lapply(v, NA-if-null))` on a JSON/YAML list, as an eval vector literal.

    Supported: all bool -> lgl; all int32 -> int; ints/floats (or out-of-int32) -> dbl; all str -> chr;
    all null -> lgl. `empty` is the vector kind of `[]` in the source dialect (None: ambiguous).
    Anything else -> literal:mixed.
    """
    if not values:
        if empty is None:
            raise Unconvertible("literal:ambiguous-empty")
        return vec(empty, [])
    nn = [v for v in values if v is not None]
    if not nn:
        return vec("lgl", [None] * len(values))
    if any(isinstance(v, (list, dict)) for v in nn):
        raise Unconvertible("literal:mixed", "nested")
    if all(isinstance(v, bool) for v in nn):
        return vec("lgl", list(values))
    if all(isinstance(v, str) for v in nn):
        return vec("chr", list(values))
    if all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in nn):
        if all(is_int32(v) for v in nn):
            return vec("int", list(values))
        return vec("dbl", [dbl_value(v) for v in values])
    raise Unconvertible("literal:mixed", f"{sorted({type(v).__name__ for v in nn})}")


def scalar_arg(v: Any) -> Any:
    """Bare scalar in a call (CONVERTERS.md, 'scalars')."""
    if isinstance(v, bool) or isinstance(v, str):
        return v
    if isinstance(v, int):
        return v if is_int32(v) else {"dbl": [float(v)]}
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v) or v.is_integer():
            return {"dbl": [dbl_value(v)]}
        return v
    raise Unconvertible("literal:scalar", repr(v))


def as_integer(values: list, ctx: str = "") -> list:
    """R as.integer() for the cases we can mirror exactly: ints, bools, integral-or-truncatable floats, NA."""
    out = []
    for v in values:
        if v is None:
            out.append(None)
        elif isinstance(v, bool):
            out.append(int(v))
        elif isinstance(v, int):
            if not is_int32(v):
                raise Unconvertible("literal:coerce", f"as.integer overflow {v} {ctx}")
            out.append(v)
        elif isinstance(v, float):
            if math.isnan(v) or math.isinf(v) or abs(v) > INT32:
                raise Unconvertible("literal:coerce", f"as.integer of {v} {ctx}")
            out.append(int(v))
        else:
            raise Unconvertible("literal:coerce", f"as.integer of {type(v).__name__} {ctx}")
    return out


class Builder:
    """Collects inputs for one case and hoists object literals out of calls."""

    def __init__(self) -> None:
        self.inputs: dict[str, dict] = {}
        self._n = 0

    def hoist(self, literal: dict, prefix: str = "_a") -> dict:
        self._n += 1
        name = f"{prefix}{self._n}"
        self.inputs[name] = literal
        return {"$ref": name}

    def fresh_name(self, prefix: str) -> str:
        self._n += 1
        return f"{prefix}{self._n}"


OBJECT_KINDS = ("IRanges", "GRanges", "Seqinfo", "DataFrame", "Factor", "IntegerList")


def literal_kind(lit: dict) -> str | None:
    if isinstance(lit, dict) and len(lit) == 1:
        k = next(iter(lit))
        if k in OBJECT_KINDS:
            return k
    return None


# ---------------------------------------------------------------- constructor trees

def _call(op: str, args: list, kwargs: dict) -> dict:
    return {"op": op, "args": args, "kwargs": kwargs}


def seqinfo_call(body: dict) -> dict:
    kw = {}
    for src, dst in (("seqnames", "seqnames"), ("seqlengths", "seqlengths"),
                     ("isCircular", "isCircular"), ("genome", "genome")):
        if src in body:
            kw[dst] = body[src]
    return _call("Seqinfo", [], kw)


def iranges_call(body: dict) -> dict:
    kw = {k: body[k] for k in ("start", "end", "width", "names") if k in body}
    return _call("IRanges", [], kw)


def constructor_tree(lit: dict, b: Builder) -> dict:
    """identity(<object literal>) -> the call that builds it (CONVERTERS.md, 'identity rewrite')."""
    kind = literal_kind(lit)
    body = lit[kind]
    if kind == "IRanges":
        tree = iranges_call(body)
    elif kind == "GRanges":
        kw: dict = {"seqnames": body["seqnames"], "ranges": {"$call": iranges_call(body)}}
        if "strand" in body:
            kw["strand"] = body["strand"]
        si = body.get("seqinfo")
        if si is not None:
            kw["seqinfo"] = si if "$ref" in si else {"$call": seqinfo_call(si)}
        tree = _call("GRanges", [], kw)
    elif kind == "Seqinfo":
        return seqinfo_call(body)
    else:
        raise Unconvertible("op-out-of-vocab:identity", f"on {kind}")
    if body.get("mcols"):
        df = b.hoist({"DataFrame": {"columns": body["mcols"]}}, "_df")
        tree = _call("mcols<-", [{"$call": tree}], {"value": df})
    return tree


def check_call(call: dict, inputs: dict) -> None:
    """Op vocabulary and ref resolution; raises Unconvertible."""
    ops = load_ops()
    for op in ops_in_call(call):
        if op not in ops:
            raise Unconvertible(f"op-out-of-vocab:{op}")

    def refs(node):
        if isinstance(node, dict):
            if "$ref" in node and len(node) == 1:
                if node["$ref"] not in inputs:
                    raise Unconvertible("unresolved-ref", node["$ref"])
            else:
                for v in node.values():
                    refs(v)
        elif isinstance(node, list):
            for v in node:
                refs(v)

    refs(call)


def prune_inputs(inputs: dict, call: dict) -> dict:
    """Keep inputs the call references, transitively, in source order."""
    need: set[str] = set()

    def walk(node):
        if isinstance(node, dict):
            if "$ref" in node and len(node) == 1:
                need.add(node["$ref"])
            else:
                for v in node.values():
                    walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)

    walk(call)
    while True:
        before = len(need)
        for nm in list(need):
            if nm in inputs:
                walk(inputs[nm])
        if len(need) == before:
            break
    return {k: v for k, v in inputs.items() if k in need}


def make_case(*, suite: str, source: str, stem: str, file: str, oid: str, tags: list[str],
              inputs: dict, call: dict) -> dict:
    check_call(call, inputs)
    return {
        "id": f"{suite}/{stem}/{oid}",
        "source": source,
        "origin": {"suite": suite, "file": file, "id": oid, "tags": tags},
        "op": call["op"],
        "inputs": prune_inputs(inputs, call),
        "call": call,
    }
