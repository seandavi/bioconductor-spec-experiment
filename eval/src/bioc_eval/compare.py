"""Comparator (spec/COMPARATOR.md): outcome per (case, port). Exact by default; relaxations R1, R2, R3 only."""

from __future__ import annotations

import math
from typing import Any

from .common import SPEC

OUTCOMES = ("pass", "mismatch", "adapter_error", "unsupported", "port_error_unexpected")

_projection_validator = None


def projection_errors(proj: Any) -> list[str]:
    """Schema errors of a port-supplied projection (a bad projection is the adapter's fault)."""
    global _projection_validator
    if _projection_validator is None:
        import json
        from jsonschema import Draft202012Validator

        _projection_validator = Draft202012Validator(json.loads((SPEC / "projection.schema.json").read_text()))
    return [e.message[:160] for e in _projection_validator.iter_errors(proj)][:3]


# ------------------------------------------------------------------ exact equality

def _is_num(x: Any) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool)


def _leaf_equal(r: Any, p: Any) -> bool:
    if r is None or p is None:
        return r is None and p is None
    if isinstance(r, bool) or isinstance(p, bool):
        return isinstance(r, bool) and isinstance(p, bool) and r == p
    if _is_num(r) and _is_num(p):
        return r == p  # -0 == 0; int 3 == dbl 3.0
    return type(r) is type(p) and r == p  # str ("NaN"/"Inf" are plain strings)


def _ulps(a: float, b: float) -> float:
    return abs(a - b) / math.ulp(max(abs(a), abs(b))) if (a != b and math.isfinite(a) and math.isfinite(b)) else 0.0


def diffs(r: Any, p: Any, path: str = "$", key: str = "") -> list[dict]:
    """All differing leaves/nodes between R and port values."""
    out: list[dict] = []
    if isinstance(r, dict) and isinstance(p, dict):
        if r.get("kind") != p.get("kind") and ("kind" in r or "kind" in p):
            return [{"path": path, "kind": "kind", "r": r.get("kind"), "port": p.get("kind")}]
        if set(r) != set(p):
            return [{"path": path, "kind": "value", "r": sorted(r), "port": sorted(p)}]
        for k in r:
            out += diffs(r[k], p[k], f"{path}.{k}", k)
        return out
    if isinstance(r, list) and isinstance(p, list):
        if len(r) != len(p):
            return [{"path": path, "kind": "length", "r": len(r), "port": len(p)}]
        scalar = all(not isinstance(x, (list, dict)) for x in r) and all(not isinstance(x, (list, dict)) for x in p)
        sub: list[dict] = []
        for i, (a, b) in enumerate(zip(r, p)):
            sub += diffs(a, b, f"{path}[{i}]", key)
        if sub and scalar and sorted(map(repr, r)) == sorted(map(repr, p)):
            return [{"path": path, "kind": "order", "r": r[:6], "port": p[:6]}]
        return sub
    if isinstance(r, (list, dict)) or isinstance(p, (list, dict)):
        return [{"path": path, "kind": "names" if key in ("names", "rownames") else "kind",
                 "r": _short(r), "port": _short(p)}]
    if _leaf_equal(r, p):
        return out
    if key == "type":
        kind = "type"
    elif key in ("names", "rownames"):
        kind = "names"
    elif _is_num(r) and _is_num(p) and (isinstance(r, float) or isinstance(p, float)) and _ulps(float(r), float(p)) <= 4:
        kind = "float_ulp"
    else:
        kind = "value"
    return [{"path": path, "kind": kind, "r": r, "port": p}]


def _short(v: Any) -> Any:
    s = repr(v)
    return s if len(s) <= 80 else s[:77] + "..."


def mismatch_detail(ds: list[dict]) -> dict:
    """First differing path + kinds over all differing leaves (float_ulp only if every difference is)."""
    kinds = {d["kind"] for d in ds}
    if "float_ulp" in kinds and kinds != {"float_ulp"}:
        kinds = (kinds - {"float_ulp"}) | {"value"}
    first = ds[0]
    return {"path": first["path"], "r": _short(first["r"]), "port": _short(first["port"]), "kinds": sorted(kinds)}


# ------------------------------------------------------------------ relaxations

def hits_query_set(r: dict, p: dict) -> list[dict]:
    """R2: query order exact, subject order within a query free. Returns [] iff equal under R2."""
    if p.get("kind") != "Hits":
        return [{"path": "$", "kind": "kind", "r": "Hits", "port": p.get("kind")}]
    ds: list[dict] = []
    for k in ("nLnode", "nRnode"):
        if r[k] != p[k]:
            ds.append({"path": f"$.{k}", "kind": "value", "r": r[k], "port": p[k]})
    pf, pt = p["from"], p["to"]
    if len(pf) != len(pt):
        return ds + [{"path": "$.to", "kind": "length", "r": len(r["to"]), "port": len(pt)}]
    if any(a > b for a, b in zip(pf, pf[1:])):
        ds.append({"path": "$.from", "kind": "order", "r": "non-decreasing", "port": "not sorted by query"})
    r_pairs = sorted(zip(r["from"], r["to"]))
    p_pairs = sorted(zip(pf, pt))
    if r_pairs != p_pairs:
        kind = "length" if len(r_pairs) != len(p_pairs) else "value"
        i = next((i for i, (a, b) in enumerate(zip(r_pairs, p_pairs)) if a != b), 0)
        ds.append({"path": "$.pairs", "kind": kind,
                   "r": r_pairs[i] if i < len(r_pairs) else len(r_pairs),
                   "port": p_pairs[i] if i < len(p_pairs) else len(p_pairs)})
    return ds


def arbitrary_in_universe(r: dict, p: dict, universe: dict) -> list[dict]:
    """R3: per query, the port's pick is one of R's select="all" hits or R's own pick (null iff none)."""
    if p.get("kind") != "vector" or p.get("type") != r["type"]:
        return [{"path": "$.type", "kind": "type" if p.get("kind") == "vector" else "kind",
                 "r": r.get("type"), "port": p.get("type", p.get("kind"))}]
    if len(p["values"]) != len(r["values"]):
        return [{"path": "$.values", "kind": "length", "r": len(r["values"]), "port": len(p["values"])}]
    if p["names"] != r["names"]:
        return [{"path": "$.names", "kind": "names", "r": r["names"], "port": p["names"]}]
    hits: dict[int, set[int]] = {}
    for f, t in zip(universe["from"], universe["to"]):
        hits.setdefault(f, set()).add(t)
    for i, v in enumerate(r["values"], start=1):
        if v is not None:  # R's own pick is always allowed (GRanges nearest select="all" can omit its own tie)
            hits.setdefault(i, set()).add(v)
    ds = []
    for i, v in enumerate(p["values"], start=1):
        allowed = hits.get(i)
        ok = (v is None) if not allowed else (v is not None and not isinstance(v, bool) and v in allowed)
        if not ok:
            ds.append({"path": f"$.values[{i - 1}]", "kind": "not_in_universe",
                       "r": sorted(allowed) if allowed else None, "port": v})
    return ds


# ------------------------------------------------------------------ outcome per (case, port)

def score(exp: dict, rec: dict | None) -> dict:
    """exp: an `ok` expected.jsonl record; rec: the port's result record (None = not run)."""
    base = {"id": exp["id"], "source": exp["source"], "op": exp["op"]}
    if rec is None:
        return {**base, "outcome": "unsupported", "rule": "not_run", "detail": {"reason": "no result record"}}
    status = rec["status"]
    r_error = exp["result"]["kind"] == "error"
    if status == "unsupported":
        return {**base, "outcome": "unsupported", "rule": "row1", "detail": {"reason": rec.get("reason", "")}}
    if status == "adapter_error":
        return {**base, "outcome": "adapter_error", "rule": "row2", "detail": rec.get("error", {})}
    if status == "port_error" and rec.get("phase") == "build":
        return {**base, "outcome": "port_error_unexpected", "rule": "row2b", "detail": rec.get("error", {})}
    if status == "port_error" and rec.get("phase") == "project":
        # the operation returned; failing to read the result back is not "the port raised" (no R1 pass)
        return {**base, "outcome": "port_error_unexpected", "rule": "row2c", "detail": rec.get("error", {})}
    if r_error and status == "port_error" and rec.get("error", {}).get("origin") != "timeout":
        return {**base, "outcome": "pass", "rule": "R1"}
    if r_error and status == "value":
        bad = projection_errors(rec["result"])
        if bad:
            return {**base, "outcome": "adapter_error", "rule": "invalid_projection", "detail": {"schema": bad}}
        return {**base, "outcome": "mismatch", "rule": "row4",
                "detail": {"path": "$", "r": "error", "port": rec["result"].get("kind"), "kinds": ["r_error_port_value"]}}
    if status == "port_error":
        return {**base, "outcome": "port_error_unexpected", "rule": "row5", "detail": rec.get("error", {})}
    p = rec["result"]
    bad = projection_errors(p)
    if bad:
        return {**base, "outcome": "adapter_error", "rule": "invalid_projection", "detail": {"schema": bad}}
    r = exp["result"]
    if "universe" in exp and r["kind"] == "vector":
        rule, ds = "R3", arbitrary_in_universe(r, p, exp["universe"])
    elif r["kind"] == "Hits" and exp.get("sorted_by_query"):
        rule, ds = "R2", hits_query_set(r, p)
    else:
        rule, ds = "exact", diffs(r, p)
    if not ds:
        return {**base, "outcome": "pass", "rule": rule}
    return {**base, "outcome": "mismatch", "rule": rule, "detail": mismatch_detail(ds)}
