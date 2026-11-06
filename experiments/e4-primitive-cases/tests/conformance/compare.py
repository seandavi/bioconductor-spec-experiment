"""Compare modes (PLAN 4.3): exact (default), ignore, oneof, multiset (rows), set."""

from __future__ import annotations

import math
from typing import Any

SKIP_KEYS = {"kind", "warnings", "message", "candidates", "type"}


def _flatten(d: dict[str, Any]) -> dict[str, Any]:
    flat: dict[str, Any] = {}
    for k, v in d.items():
        if k in SKIP_KEYS:
            continue
        if k in ("mcols", "seqinfo"):
            flat.update({f"{k}.{c}": cv for c, cv in v.items()})
        else:
            flat[k] = v
    return flat


def _eq(a: Any, b: Any) -> bool:
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_eq(x, y) for x, y in zip(a, b, strict=True))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_eq(a[k], b[k]) for k in a)
    if a is None or b is None:
        return a is b
    if isinstance(a, bool) != isinstance(b, bool):
        return False
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return a == b or math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)
    return bool(a == b)


def _canon(v: Any) -> Any:
    return sorted((_canon(e) for e in v), key=repr) if isinstance(v, list) else v


def compare(spec: dict[str, Any], golden: dict[str, Any], actual: dict[str, Any]) -> list[str]:
    """Return mismatch descriptions (empty = pass)."""
    cmp = spec.get("compare") or {}
    modes: dict[str, str] = cmp.get("fields") or {}
    problems: list[str] = []
    if golden["kind"] != actual["kind"]:
        return [f"kind: golden {golden['kind']} != actual {actual['kind']}"]
    g, a = _flatten(golden), _flatten(actual)
    if g.keys() != a.keys():
        problems.append(f"fields differ: golden {sorted(g)} vs actual {sorted(a)}")
    shared = [k for k in g if k in a and modes.get(k) != "ignore"]
    if cmp.get("rows") == "multiset":
        row_fields = [k for k in shared if isinstance(g[k], list) and not k.startswith("seqinfo.")]
        gr = sorted(map(repr, zip(*(g[k] for k in row_fields), strict=False)))
        ar = sorted(map(repr, zip(*(a[k] for k in row_fields), strict=False)))
        if gr != ar:
            problems.append(f"rows (multiset over {row_fields}): golden {gr} != actual {ar}")
        shared = [k for k in shared if k not in row_fields]
    for k in shared:
        mode = modes.get(k, "exact")
        if mode == "exact":
            ok = _eq(g[k], a[k])
        elif mode == "set":
            ok = _eq(_canon(g[k]), _canon(a[k]))
        elif mode == "oneof":
            cands = golden["candidates"]
            ok = len(a[k]) == len(cands) and all(
                any(_eq(v, c) for c in cs) for v, cs in zip(a[k], cands, strict=True)
            )
        else:
            raise ValueError(f"unknown compare mode {mode!r}")
        if not ok:
            problems.append(f"{k} [{mode}]: golden {g[k]!r} != actual {a[k]!r}")
    return problems
