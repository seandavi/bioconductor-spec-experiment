"""E1 (differential suite, YAML in R vocabulary) -> eval cases. Mirrors E1's own R decoder (r/run_cases.R)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from ..common import WORKTREES
from .common import (
    Builder, Unconvertible, as_integer, constructor_tree, dbl_value, is_int32, literal_kind,
    make_case, scalar_arg, untyped_vector, vec,
)

ROOT = WORKTREES / "e1-differential/experiments/e1-differential"
CASES_DIR = ROOT / "tests/differential/cases"
SUITE_PREFIX = "tests/differential/cases"

# the four `raw` cases, mapped by hand (CONVERTERS.md); inputs: none
RAW = {
    "rt-edge-ir-default-empty": {"op": "IRanges", "args": [], "kwargs": {}},
    "rt-edge-gr-default-empty": {"op": "GRanges", "args": [], "kwargs": {}},
    "rt-meta-seqinfo-numeric-names": {"op": "Seqinfo", "args": [{"int": [1, 2]}], "kwargs": {}},
    "rt-meta-gr-recycle-range-names": {
        "op": "GRanges",
        "args": [{"chr": ["a", "b"]},
                 {"$call": {"op": "IRanges", "args": [1], "kwargs": {"width": 2, "names": "only"}}}],
        "kwargs": {},
    },
}


def _typed(kind: str, v: Any) -> dict:
    """`{int|num|chr|lgl: [...]}` -> typed vector (R: unlist, then as.<type>)."""
    vals = v if isinstance(v, list) else [v]
    if kind == "int":
        return vec("int", as_integer(vals))
    if kind == "num":
        out = []
        for e in vals:
            if e is not None and (isinstance(e, str) or isinstance(e, (list, dict))):
                raise Unconvertible("literal:coerce", "as.numeric of non-number")
            out.append(dbl_value(e if not isinstance(e, bool) else float(e)))
        return vec("dbl", out)
    if kind == "chr":
        if any(e is not None and not isinstance(e, str) for e in vals):
            raise Unconvertible("literal:coerce", "as.character of non-string")
        return vec("chr", list(vals))
    if any(e is not None and not isinstance(e, bool) for e in vals):
        raise Unconvertible("literal:coerce", "as.logical of non-bool")
    return vec("lgl", list(vals))


def vector_literal(v: Any) -> dict:
    """decode() result that must be an atomic vector."""
    if v is None:
        return vec("lgl", [None])
    if isinstance(v, dict):
        k = next(iter(v))
        if len(v) > 1:
            # E1's R decoder reads only the first key; the extra keys are a source typo (misplaced brace)
            raise Unconvertible("literal:source-multi-key", ",".join(v))
        if k in ("int", "num", "chr", "lgl"):
            return _typed(k, v[k])
        raise Unconvertible("literal:not-a-vector", k)
    if isinstance(v, list):
        return untyped_vector(v, empty="int")
    sc = scalar_arg(v)
    if isinstance(sc, dict):
        return sc
    kind = "lgl" if isinstance(sc, bool) else "chr" if isinstance(sc, str) else "int" if isinstance(sc, int) else "dbl"
    return vec(kind, [sc])


def columns(cols: Any) -> list:
    if not isinstance(cols, dict) or not cols:
        raise Unconvertible("literal:empty-mcols" if not cols else "literal:columns")
    out = []
    for name, v in cols.items():
        out.append([str(name), vector_literal(v)])
    return out


def seqinfo_body(l: dict) -> dict:
    if l.get("seqnames") is None:
        raise Unconvertible("literal:seqinfo-no-seqnames")
    body = {"seqnames": vector_literal(l["seqnames"])}
    if l.get("seqlengths") is not None:
        sl = vector_literal(l["seqlengths"])
        kind = next(k for k in sl if k != "names")
        if kind != "int":
            if kind not in ("dbl", "lgl"):
                raise Unconvertible("literal:coerce", "as.integer of seqlengths")
            sl = vec("int", as_integer(sl[kind]))
        body["seqlengths"] = sl
    if l.get("is_circular") is not None:
        body["isCircular"] = vector_literal(l["is_circular"])
    if l.get("genome") is not None:
        body["genome"] = vector_literal(l["genome"])
    return body


def ranges_body(l: dict, gr: bool) -> dict:
    allowed = {"start", "end", "width", "names", "mcols"} | ({"seqnames", "strand", "seqinfo"} if gr else set())
    extra = set(l) - allowed
    if extra:
        raise Unconvertible("literal:unknown-field", ",".join(sorted(extra)))
    body: dict = {}
    if gr:
        if l.get("seqnames") is None:
            raise Unconvertible("literal:granges-no-seqnames")
        body["seqnames"] = vector_literal(l["seqnames"])
    for k in ("start", "end", "width", "names"):
        if l.get(k) is not None and not (isinstance(l[k], dict) and list(l[k]) == ["NULL"]):
            body[k] = vector_literal(l[k])
    if gr:
        if l.get("strand") is not None:
            body["strand"] = vector_literal(l["strand"])
        if l.get("seqinfo") is not None:
            if not isinstance(l["seqinfo"], dict):
                raise Unconvertible("literal:seqinfo-ref")
            body["seqinfo"] = seqinfo_body(l["seqinfo"])
    if l.get("mcols") is not None:
        body["mcols"] = columns(l["mcols"])
    return body


def object_literal(v: dict) -> dict | None:
    """IRanges / GRanges / Seqinfo / DataFrame literal -> eval literal, else None."""
    if len(v) != 1:
        return None
    k = next(iter(v))
    body = v[k]
    if k == "IRanges":
        return {"IRanges": ranges_body(body, False)}
    if k == "GRanges":
        return {"GRanges": ranges_body(body, True)}
    if k == "Seqinfo":
        return {"Seqinfo": seqinfo_body(body)}
    if k == "DataFrame":
        return {"DataFrame": {"columns": columns(body)}}
    return None


def input_literal(v: Any) -> dict:
    if isinstance(v, dict):
        if list(v) == ["NULL"]:
            raise Unconvertible("literal:null-input")
        if (lit := object_literal(v)) is not None:
            return lit
    return vector_literal(v)


def convert_arg(v: Any, b: Builder) -> Any:
    if isinstance(v, str):
        if v.startswith("$$"):
            return v[1:]
        if v.startswith("$"):
            return {"$ref": v[1:]}
        return v
    if v is None:
        return vec("lgl", [None])
    if isinstance(v, dict):
        if list(v) == ["call"]:
            return {"$call": convert_call(v["call"], b)}
        if list(v) == ["NULL"]:
            return None
        if (lit := object_literal(v)) is not None:
            return b.hoist(lit)
        return vector_literal(v)
    if isinstance(v, list):
        return untyped_vector(v, empty="int")
    return scalar_arg(v)


def convert_call(c: dict, b: Builder) -> dict:
    extra = set(c) - {"fn", "args", "kwargs"}
    if extra:
        raise Unconvertible("call:unknown-field", ",".join(sorted(extra)))
    args = [convert_arg(a, b) for a in c.get("args") or []]
    kwargs = {k: convert_arg(v, b) for k, v in (c.get("kwargs") or {}).items()}
    return {"op": c["fn"], "args": args, "kwargs": kwargs}


def convert_case(spec: dict, stem: str, file: str) -> dict:
    cid = spec["id"]
    tags = [str(t) for t in spec.get("tags") or []]
    source = "e1-redteam" if stem.startswith("redteam_") or "redteam" in tags else "e1"
    b = Builder()
    if "raw" in spec:
        if cid not in RAW:
            raise Unconvertible("raw-r", spec["raw"].get("r", ""))
        call = RAW[cid]
        return make_case(suite="e1", source=source, stem=stem, file=file, oid=cid, tags=tags, inputs={}, call=call)
    inputs = {name: input_literal(v) for name, v in (spec.get("inputs") or {}).items()}
    b.inputs = dict(inputs)
    c = spec["call"]
    if c["fn"] == "identity":
        args = c.get("args") or []
        if (len(args) == 1 and not c.get("kwargs") and isinstance(args[0], str) and args[0].startswith("$")
                and args[0][1:] in inputs and literal_kind(inputs[args[0][1:]]) in ("IRanges", "GRanges", "Seqinfo")):
            call = constructor_tree(inputs[args[0][1:]], b)
        else:
            raise Unconvertible("op-out-of-vocab:identity")
    else:
        call = convert_call(c, b)
    return make_case(suite="e1", source=source, stem=stem, file=file, oid=cid, tags=tags, inputs=b.inputs, call=call)


def convert() -> tuple[list[dict], list[dict], dict]:
    cases, unconv = [], []
    read = 0
    for path in sorted(CASES_DIR.glob("*.yaml")):
        doc = yaml.safe_load(path.read_text()) or {}
        file = f"{SUITE_PREFIX}/{path.name}"
        for spec in doc.get("cases", []):
            read += 1
            try:
                cases.append(convert_case(spec, path.stem, file))
            except Unconvertible as e:
                unconv.append({"source": "e1", "origin": {"suite": "e1", "file": file, "id": spec["id"],
                                                          "tags": [str(t) for t in spec.get("tags") or []]},
                               "reason": e.reason, "detail": e.detail[:200]})
    return cases, unconv, {"read": read}
