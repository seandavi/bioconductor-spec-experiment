"""Print clusters and their examined cases (R value, port value, inputs) for manual adjudication.

Usage: python adjudication/show.py PORT [SIGNATURE_SUBSTRING] [K]
Examined cases are `picks(ids, K)`: K ids spread evenly over the cluster (all when n <= K).
"""
from __future__ import annotations

import collections
import json
import sys

from load import load
from cluster import OUT


def picks(ids: list[str], k: int = 3) -> list[str]:
    ids = sorted(ids)
    return ids if len(ids) <= k else [ids[i * len(ids) // k] for i in range(k)]


def _compact(c: dict) -> str:
    """Inputs summarised to their fields (lengths for long vectors), then the call."""
    def lit(v):
        if isinstance(v, dict) and v.get("kind") in ("IRanges", "GRanges"):
            keys = ("seqnames", "start", "width", "strand", "names")
            body = {k: v[k] for k in keys if k in v and v[k] is not None}
            if v.get("mcols"):
                body["mcols"] = [n for n, _ in v["mcols"]]
            if "seqinfo" in v:
                si = v["seqinfo"]
                body["si"] = [si["seqnames"], si["seqlengths"], si["isCircular"]]
            return f"{v['kind']}{json.dumps(body)}"
        return json.dumps(v)
    ins = "; ".join(f"{k}={lit(v)}" for k, v in c["inputs"].items())
    return f"{ins}\n     CALL {json.dumps(c['call'])}"


def _short(o, n: int) -> str:
    s = json.dumps(o)
    return s if len(s) < n else s[:n] + "..."


def main() -> None:
    port = sys.argv[1]
    sub = sys.argv[2] if len(sys.argv) > 2 else None
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    resolved, expected, raw, scored = load()
    by = collections.defaultdict(list)
    for r in map(json.loads, OUT.open()):
        if r["port"] == port:
            by[r["signature"]].append(r)
    for sig, rs in sorted(by.items(), key=lambda x: (-len(x[1]), x[0])):
        if sub and sub not in sig:
            continue
        ops = collections.Counter(r["op"] for r in rs)
        print(f"\n######## {len(rs)} {sig}  ops={dict(ops)}")
        if sub is None:
            continue
        for cid in picks([r["id"] for r in rs], k):
            c, e, p = resolved[cid], expected[cid], raw[port][cid]
            print("--", cid, scored[port][cid].get("detail"))
            print("  IN", _short(_compact(c), 1400))
            print("  R ", _short(e["result"], 800))
            print("  P ", _short(p.get("result") or p.get("error"), 800))


if __name__ == "__main__":
    main()
