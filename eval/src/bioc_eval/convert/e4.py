"""E4 (primitive-cases suite, YAML) -> eval cases. Mirrors E4's own R builder (runners/r/build.R)."""

from __future__ import annotations

from typing import Any

import yaml

from ..common import WORKTREES
from .common import (
    Builder, Unconvertible, as_integer, constructor_tree, dbl_value, literal_kind, make_case, scalar_arg,
    untyped_vector, vec,
)

ROOT = WORKTREES / "e4-primitive-cases/experiments/e4-primitive-cases"
CASES_DIR = ROOT / "cases"


def untyped(v: Any) -> dict:
    """E4 `vec(x)` with no type: R yaml simplification, then unlist. `[]` is logical(0)."""
    if isinstance(v, list):
        return untyped_vector(v, empty="lgl")
    sc = scalar_arg(v)
    if isinstance(sc, dict):
        return sc
    return vec("lgl" if isinstance(sc, bool) else "chr" if isinstance(sc, str) else "int" if isinstance(sc, int) else "dbl", [sc])


def typed(v: Any, kind: str) -> dict:
    """E4 `vec(x, type)`: unlist then as.<type>."""
    vals = v if isinstance(v, list) else [v]
    if kind == "int":
        return vec("int", as_integer(vals))
    if kind == "num":
        if any(isinstance(e, str) for e in vals):
            raise Unconvertible("literal:coerce", "as.numeric of string")
        return vec("dbl", [dbl_value(float(e) if isinstance(e, bool) else e) for e in vals])
    if kind == "chr":
        if any(e is not None and not isinstance(e, str) for e in vals):
            raise Unconvertible("literal:coerce", "as.character of non-string")
        return vec("chr", list(vals))
    if any(e is not None and not isinstance(e, bool) for e in vals):
        raise Unconvertible("literal:coerce", "as.logical of non-bool")
    return vec("lgl", list(vals))


def mcols(mc: Any) -> list:
    if not isinstance(mc, dict) or not mc:
        raise Unconvertible("literal:empty-mcols")
    return [[str(k), untyped(v)] for k, v in mc.items()]


def seqinfo_body(s: dict) -> dict:
    body = {"seqnames": typed(s.get("seqnames"), "chr")}
    if s.get("seqlengths") is not None:
        body["seqlengths"] = typed(s["seqlengths"], "int")
    if s.get("isCircular") is not None:
        body["isCircular"] = typed(s["isCircular"], "lgl")
    if s.get("genome") is not None:
        body["genome"] = typed(s["genome"], "chr")
    return body


def ranges_body(s: dict, with_mcols: bool = True) -> dict:
    body = {}
    for k in ("start", "end", "width"):
        if s.get(k) is not None:
            body[k] = typed(s[k], "int")
    if s.get("names") is not None:
        body["names"] = typed(s["names"], "chr")
    if with_mcols and s.get("mcols") is not None:
        body["mcols"] = mcols(s["mcols"])
    return body


def input_literal(spec: dict, inputs: dict) -> dict:
    if len(spec) != 1:
        raise Unconvertible("input:shape")
    kind = next(iter(spec))
    s = spec[kind]
    if kind == "IRanges":
        return {"IRanges": ranges_body(s)}
    if kind == "Seqinfo":
        return {"Seqinfo": seqinfo_body(s)}
    if kind == "GRanges":
        body = {"seqnames": typed(s.get("seqnames"), "chr")}
        body.update(ranges_body(s, with_mcols=False))
        if s.get("strand") is not None:
            body["strand"] = typed(s["strand"], "chr")
        si = s.get("seqinfo")
        if isinstance(si, str):
            if si not in inputs:
                raise Unconvertible("unresolved-ref", si)
            body["seqinfo"] = {"$ref": si}
        elif si is not None:
            body["seqinfo"] = seqinfo_body(si)
        if s.get("mcols") is not None:
            body["mcols"] = mcols(s["mcols"])
        return {"GRanges": body}
    if kind in ("int", "num", "chr", "lgl"):
        return typed(s, kind)
    raise Unconvertible(f"input-kind:{kind}")


def kwarg(v: Any) -> Any:
    if isinstance(v, dict):
        if list(v) == ["ref"]:
            return {"$ref": v["ref"]}
        raise Unconvertible("kwarg:dict")
    if v is None:
        return None
    if isinstance(v, list):
        return untyped(v)
    return scalar_arg(v)


def convert_case(spec: dict, domain: str, stem: str, file: str) -> dict:
    tags = [domain, str(spec.get("clause"))]
    b = Builder()
    for name, ispec in spec["inputs"].items():
        b.inputs[name] = input_literal(ispec, b.inputs)
    c = spec["call"]
    if c["op"] == "identity":
        names = c.get("args") or []
        if (len(names) == 1 and not c.get("kwargs") and names[0] in b.inputs
                and literal_kind(b.inputs[names[0]]) in ("IRanges", "GRanges", "Seqinfo")):
            call = constructor_tree(b.inputs[names[0]], b)
        else:
            raise Unconvertible("op-out-of-vocab:identity")
    else:
        call = {"op": c["op"], "args": [{"$ref": n} for n in c.get("args") or []],
                "kwargs": {k: kwarg(v) for k, v in (c.get("kwargs") or {}).items()}}
    return make_case(suite="e4", source="e4", stem=f"{domain}-{stem}", file=file, oid=spec["id"], tags=tags,
                     inputs=b.inputs, call=call)


def convert() -> tuple[list[dict], list[dict], dict]:
    cases, unconv, read = [], [], 0
    for path in sorted(CASES_DIR.glob("*/*.yaml")):
        domain = path.parent.name
        file = f"cases/{domain}/{path.name}"
        doc = yaml.safe_load(path.read_text()) or {}
        for spec in doc.get("cases", []):
            read += 1
            try:
                cases.append(convert_case(spec, domain, path.stem, file))
            except Unconvertible as e:
                unconv.append({"source": "e4", "origin": {"suite": "e4", "file": file, "id": spec["id"],
                                                          "tags": [domain, str(spec.get("clause"))]},
                               "reason": e.reason, "detail": e.detail[:200]})
    return cases, unconv, {"read": read}
