"""E2 (conformance suite, JSON with $-tags) -> eval cases. Mirrors E2's own R decoder (conformance/r/codec.R)."""

from __future__ import annotations

import json
from typing import Any

from ..common import WORKTREES
from .common import (
    Builder, Unconvertible, as_integer, dbl_value, is_int32, make_case, scalar_arg, untyped_vector, vec,
)

ROOT = WORKTREES / "e2-conformance/experiments/e2-conformance"
CASES_DIR = ROOT / "conformance/cases"
SUITE_PREFIX = "conformance/cases"
TAGS = {"$int": "int", "$dbl": "dbl", "$chr": "chr", "$lgl": "lgl"}


class Decoder:
    """Decodes one file's fixtures and case arguments into eval literals."""

    def __init__(self, fixtures: dict):
        self.fixtures = fixtures
        self.lits: dict[str, dict] = {}

    # ---- vectors
    def tagged(self, tree: dict) -> dict:
        tag = next(t for t in tree if t in TAGS)
        kind, items = TAGS[tag], tree[tag]
        if kind == "int":
            vals = as_integer(items)
        elif kind == "dbl":
            vals = []
            for e in items:
                if e is not None and isinstance(e, str) and e not in ("NaN", "Inf", "-Inf"):
                    raise Unconvertible("literal:coerce", f"as.numeric({e!r})")
                vals.append(e if isinstance(e, str) else dbl_value(e))
        elif kind == "chr":
            if any(e is not None and not isinstance(e, str) for e in items):
                raise Unconvertible("literal:coerce", "$chr with non-string")
            vals = list(items)
        else:
            if any(e is not None and not isinstance(e, bool) for e in items):
                raise Unconvertible("literal:coerce", "$lgl with non-bool")
            vals = list(items)
        names = tree.get("names")
        if names is not None and any(n is not None and not isinstance(n, str) for n in names):
            raise Unconvertible("literal:coerce", "names")
        return vec(kind, vals, list(names) if names is not None else None)

    def shorthand(self, v: Any) -> dict:
        items = v if isinstance(v, list) else [v]
        if not items:
            raise Unconvertible("literal:ambiguous-empty")
        if any(isinstance(e, (list, dict)) for e in items):
            raise Unconvertible("literal:mixed", "nested")
        nn = [e for e in items if e is not None]
        if not nn:
            return vec("lgl", [None] * len(items))
        if all(isinstance(e, str) for e in nn):
            return vec("chr", list(items))
        if all(isinstance(e, bool) for e in nn):
            return vec("lgl", list(items))
        if all(isinstance(e, (int, float)) and not isinstance(e, bool) for e in nn):
            if any(isinstance(e, float) for e in nn) or not all(is_int32(e) for e in nn):
                return vec("dbl", [dbl_value(e) for e in items])
            return vec("int", list(items))
        raise Unconvertible("literal:mixed", "shorthand kinds")

    def vector(self, tree: Any) -> dict:
        if isinstance(tree, dict):
            return self.tagged(tree)
        return self.shorthand(tree)

    def maybe_vector(self, tree: Any) -> dict | None:
        if tree is None:
            return None
        if isinstance(tree, dict) and "$ref" in tree:
            lit = self.fixture(tree["$ref"])
            if len(lit) == 1 and next(iter(lit)) in ("int", "dbl", "chr", "lgl"):
                return lit
            raise Unconvertible("fixture-ref:vector-field", tree["$ref"])
        return self.vector(tree)

    # ---- objects
    def columns(self, tree: Any) -> list:
        """mcols DFrame (inline or $ref) -> [[name, literal]]."""
        if isinstance(tree, dict) and "$ref" in tree:
            lit = self.fixture(tree["$ref"])
            if "DataFrame" not in lit:
                raise Unconvertible("fixture-ref:mcols", tree["$ref"])
            return lit["DataFrame"]["columns"]
        if not (isinstance(tree, dict) and tree.get("$type") == "DFrame"):
            raise Unconvertible("literal:mcols-not-dframe")
        return [[name, self.value_literal(col)] for name, col in (tree.get("columns") or {}).items()]

    def ranges_fields(self, tree: dict) -> dict:
        body: dict = {}
        for k in ("start", "end", "width"):
            if (v := self.maybe_vector(tree.get(k))) is not None:
                body[k] = v
        if (v := self.maybe_vector(tree.get("names"))) is not None:
            body["names"] = v
        if tree.get("mcols") is not None:
            cols = self.columns(tree["mcols"])
            if cols:
                body["mcols"] = cols
        return body

    def seqinfo_body(self, tree: Any) -> dict:
        if isinstance(tree, dict) and "$ref" in tree:
            lit = self.fixture(tree["$ref"])
            if "Seqinfo" not in lit:
                raise Unconvertible("fixture-ref:seqinfo", tree["$ref"])
            return lit["Seqinfo"]
        if not (isinstance(tree, dict) and tree.get("$type") == "Seqinfo"):
            raise Unconvertible("literal:seqinfo")
        body = {}
        for src in ("seqnames", "seqlengths", "isCircular", "genome"):
            if (v := self.maybe_vector(tree.get(src))) is not None:
                body[src] = v
        return body

    def value_literal(self, tree: Any) -> dict:
        """Any decodable value -> eval literal (vector or object)."""
        if tree is None:
            raise Unconvertible("literal:null-value")
        if not isinstance(tree, dict):
            return self.shorthand(tree)
        if "$ref" in tree:
            return {"$ref": tree["$ref"]}
        if any(t in tree for t in TAGS):
            return self.tagged(tree)
        t = tree.get("$type")
        if t == "IRanges":
            return {"IRanges": self.ranges_fields(tree)}
        if t == "GRanges":
            sn = self.maybe_vector(tree.get("seqnames"))
            if sn is None or "chr" not in sn:
                raise Unconvertible("literal:granges-seqnames")
            body = {"seqnames": sn}
            body.update(self.ranges_fields(tree))
            # E2's decode made explicit: strand defaults to rep("*", n); no seqinfo -> Seqinfo(unique(seqnames))
            body["strand"] = self.maybe_vector(tree.get("strand")) or vec("chr", ["*"] * len(sn["chr"]))
            t_si = tree.get("seqinfo")
            if t_si is None:
                body["seqinfo"] = {"seqnames": vec("chr", list(dict.fromkeys(sn["chr"])))}
            elif isinstance(t_si, dict) and "$ref" in t_si:
                self.seqinfo_body(t_si)  # validates the fixture is a Seqinfo
                body["seqinfo"] = {"$ref": t_si["$ref"]}
            else:
                body["seqinfo"] = self.seqinfo_body(t_si)
            order = ["seqnames", "start", "end", "width", "names", "strand", "seqinfo", "mcols"]
            return {"GRanges": {k: body[k] for k in order if k in body}}
        if t == "Seqinfo":
            return {"Seqinfo": self.seqinfo_body(tree)}
        if t == "DFrame":
            cols = [[name, self.value_literal(col)] for name, col in (tree.get("columns") or {}).items()]
            if tree.get("nrow") is not None:
                n = tree["nrow"]
            else:
                n = literal_length(cols[0][1]) if cols else 0
            return {"DataFrame": {"columns": cols, "nrow": n}}
        if t == "Factor":
            lv, cd = self.vector(tree["levels"]), self.vector(tree["codes"])
            if "chr" not in lv or not ("int" in cd or "lgl" in cd):
                raise Unconvertible("literal:factor")
            codes = as_integer(next(v for k, v in cd.items() if k != "names"))
            return {"Factor": {"levels": lv["chr"], "codes": codes}}
        if t == "IntegerList":
            els = [as_integer(list(e)) for e in tree["elements"]]
            body = {"elements": els}
            if tree.get("names") is not None:
                nm = self.vector(tree["names"])
                body["names"] = nm["chr"]
            return {"IntegerList": body}
        if t in ("Hits", "Rle", "RleList"):
            raise Unconvertible(f"input-kind:{t}")
        raise Unconvertible("literal:unknown-type", str(t))

    def fixture(self, name: str) -> dict:
        if name not in self.lits:
            if name not in self.fixtures:
                raise Unconvertible("unresolved-ref", name)
            self.lits[name] = self.value_literal(self.fixtures[name])
        return self.lits[name]

    # ---- call arguments
    def arg(self, v: Any, b: Builder) -> Any:
        if v is None:
            return None
        if isinstance(v, dict):
            if "$ref" in v:
                self.fixture(v["$ref"])
                return {"$ref": v["$ref"]}
            if any(t in v for t in TAGS):
                return self.tagged(v)
            lit = self.value_literal(v)
            return b.hoist(lit)
        if isinstance(v, list):
            return self.shorthand(v)
        return scalar_arg(v)


def literal_length(lit: dict) -> int:
    if "Factor" in lit:
        return len(lit["Factor"]["codes"])
    if "IntegerList" in lit:
        return len(lit["IntegerList"]["elements"])
    return len(next(v for k, v in lit.items() if k != "names"))


def refs_in(node: Any) -> set[str]:
    out: set[str] = set()
    if isinstance(node, dict):
        if "$ref" in node:
            out.add(node["$ref"])
        for v in node.values():
            out |= refs_in(v)
    elif isinstance(node, list):
        for v in node:
            out |= refs_in(v)
    return out


def convert_file(path) -> tuple[list[dict], list[dict], int]:
    doc = json.loads(path.read_text())
    group = doc["group"]
    file = f"{SUITE_PREFIX}/{path.name}"
    cases, unconv = [], []
    for cs in doc["cases"]:
        tags = [group, cs["clause"]]
        try:
            dec = Decoder(doc["fixtures"])
            b = Builder()
            # fixtures first (source order), restricted to what the case references (transitively)
            need = refs_in([cs.get("args"), cs.get("kwargs")])
            order = list(doc["fixtures"])
            closure = set()
            stack = list(need)
            while stack:
                nm = stack.pop()
                if nm in closure:
                    continue
                if nm not in doc["fixtures"]:
                    raise Unconvertible("unresolved-ref", nm)
                closure.add(nm)
                stack.extend(refs_in(doc["fixtures"][nm]))
            for nm in order:
                if nm in closure:
                    b.inputs[nm] = dec.fixture(nm)
            args = [dec.arg(a, b) for a in cs.get("args") or []]
            kwargs = {k: dec.arg(v, b) for k, v in (cs.get("kwargs") or {}).items()}
            call = {"op": cs["op"], "args": args, "kwargs": kwargs}
            cases.append(make_case(suite="e2", source="e2", stem=path.stem, file=file, oid=cs["id"], tags=tags,
                                   inputs=b.inputs, call=call))
        except Unconvertible as e:
            unconv.append({"source": "e2", "origin": {"suite": "e2", "file": file, "id": cs["id"], "tags": tags},
                           "reason": e.reason, "detail": e.detail[:200]})
    return cases, unconv, len(doc["cases"])


def convert() -> tuple[list[dict], list[dict], dict]:
    cases, unconv, read = [], [], 0
    for path in sorted(CASES_DIR.glob("*.json")):
        c, u, n = convert_file(path)
        cases += c
        unconv += u
        read += n
    return cases, unconv, {"read": read}
