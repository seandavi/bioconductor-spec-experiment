"""R self-consistency probes for coverage() with NA weights.

Coverage at position p is the sum of the weights of the ranges covering p. Two checks per case
whose `weight` holds NA:

1. `local_na`: R's value at every position equals that definition with NA propagating only
   from the ranges that carry it (NA iff an NA-weighted range covers p). Positions are compared
   over R's own result length, so R's default width is not re-derived. `False` = R puts NA (or a
   value) where no covering range justifies it.
2. `order_dependent`: R is re-run with the ranges (and their per-range weight/shift) reversed;
   a sum cannot depend on order, so a changed result is R being inconsistent with itself.

Writes probe_coverage_na.jsonl: {id, local_na, order_dependent, r, r_reversed}.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

from load import EVAL
from rprobe import literal_cases, run

OUT = Path(__file__).resolve().parent / "probe_coverage_na.jsonl"
VEC = ("int", "dbl", "chr", "lgl")


def _rev_vec(v: dict) -> dict:
    v = dict(v)
    for k in VEC:
        if k in v:
            v[k] = v[k][::-1]
    if v.get("names") is not None:
        v["names"] = v["names"][::-1]
    return v


def _len(v) -> int | None:
    return next((len(v[k]) for k in VEC if isinstance(v, dict) and k in v), None)


def reverse_ranges(case: dict) -> dict | None:
    c = copy.deepcopy(case)
    x = c["call"]["args"][0]
    if not (isinstance(x, dict) and "$ref" in x):
        return None
    lit = c["inputs"][x["$ref"]]
    body = lit.get("IRanges") or lit.get("GRanges")
    if body is None:
        return None
    n = next(_len(body[k]) for k in ("start", "end", "width") if k in body)
    for k in ("seqnames", "start", "end", "width", "strand", "names"):
        if k in body and _len(body[k]) == n:
            body[k] = _rev_vec(body[k])
    for col in body.get("mcols", []):
        if _len(col[1]) == n:
            col[1] = _rev_vec(col[1])
    for k in ("weight", "shift"):
        v = c["call"]["kwargs"].get(k)
        if isinstance(v, dict) and _len(v) == n and n > 1:
            c["call"]["kwargs"][k] = _rev_vec(v)
    c["id"] += "#reversed"
    return c


def _recycle(v: dict | None, n: int, default):
    if v is None:
        return [default] * n
    if not isinstance(v, dict):
        return [v] * n
    vals = next(v[k] for k in VEC if k in v)
    return [vals[i % len(vals)] for i in range(n)]


def _expand(rle: dict) -> list:
    return [v for v, k in zip(rle["values"], rle["lengths"]) for _ in range(k)]


def local_na(case: dict, r: dict) -> bool | None:
    """R's coverage equals the local-NA definition at every position of R's result (None: not checkable)."""
    if r["kind"] not in ("Rle", "RleList"):
        return None
    x = case["inputs"][case["call"]["args"][0]["$ref"]]
    body = x.get("IRanges") or x.get("GRanges")
    kw = case["call"]["kwargs"]
    if any(k not in ("weight", "shift", "width", "method") for k in kw) or len(case["call"]["args"]) > 1:
        return None
    n = next(_len(body[k]) for k in ("start", "end", "width") if k in body)
    start = _recycle(body.get("start"), n, None)
    width = _recycle(body.get("width"), n, None)
    end = _recycle(body.get("end"), n, None)
    if None in start:
        start = [e - w + 1 for e, w in zip(end, width)]
    if None in width:
        width = [e - s + 1 for s, e in zip(start, end)]
    if isinstance(kw.get("weight"), str):
        return None
    weight = _recycle(kw.get("weight"), n, 1)
    shift = _recycle(kw.get("shift"), n, 0)
    if not all(isinstance(v, int) for v in shift):
        return None
    if r["kind"] == "Rle":
        groups = {None: (_expand(r), list(range(n)))}
    else:
        seq = _recycle(body["seqnames"], n, None)
        groups = {nm: (_expand(el), [i for i in range(n) if seq[i] == nm])
                  for nm, el in zip(r["names"], r["elements"])}
    for got, idx in groups.values():
        for p in range(1, len(got) + 1):
            cov = [i for i in idx if width[i] > 0 and start[i] + shift[i] <= p <= start[i] + shift[i] + width[i] - 1]
            want = None if any(weight[i] is None for i in cov) else sum(weight[i] for i in cov)
            if (want is None) != (got[p - 1] is None) or (want is not None and want != got[p - 1]):
                return False
    return True


def main() -> None:
    expected = {r["id"]: r for r in map(json.loads, (EVAL / "build" / "expected.jsonl").open())}
    cases = literal_cases()
    todo = []
    for cid, c in cases.items():
        w = c["call"]["kwargs"].get("weight")
        if c["op"] == "coverage" and isinstance(w, dict) and any(
                v is None for k in VEC for v in w.get(k, [])):
            rc = reverse_ranges(c)
            if rc is not None:
                todo.append(rc)
    got = run(todo)
    n_dep = n_nonlocal = 0
    with OUT.open("w") as fh:
        for rc in todo:
            cid = rc["id"].removesuffix("#reversed")
            r0, r1 = expected[cid]["result"], got[rc["id"]].get("result")
            dep = r0 != r1
            loc = local_na(cases[cid], r0)
            n_dep += dep
            n_nonlocal += loc is False
            fh.write(json.dumps({"id": cid, "local_na": loc, "order_dependent": dep, "r": r0, "r_reversed": r1}) + "\n")
    print(f"{len(todo)} NA-weight coverage cases: R not local-NA in {n_nonlocal}; "
          f"R changes under range reversal in {n_dep}")


if __name__ == "__main__":
    main()
