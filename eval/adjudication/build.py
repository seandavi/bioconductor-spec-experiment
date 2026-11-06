"""Apply verdicts.RULES to signatures.jsonl -> clusters.jsonl + summary.json.

Usage: python adjudication/build.py [--check]   (--check: list clusters no rule matches, exit 1 if any)
"""
from __future__ import annotations

import collections
import json
import re
import sys
from pathlib import Path

from cluster import OUT as SIGNATURES
from load import PORTS, load
from show import picks
from verdicts import RULES

HERE = Path(__file__).resolve().parent
CLASSES = ("port_bug", "declared_divergence", "r_accident", "adapter_fault", "harness_fault")
EXCLUDED = ("r_accident", "declared_divergence", "adapter_fault", "harness_fault")


def verdict(sig: str) -> tuple[int, tuple[str, str, str, str]] | None:
    for i, rule in enumerate(RULES):
        if re.search(rule[0], sig):
            return i, rule
    return None


def main() -> None:
    rows = [json.loads(line) for line in SIGNATURES.open()]
    by: dict[str, list[dict]] = collections.defaultdict(list)
    for r in rows:
        by[r["signature"]].append(r)
    missing = sorted((len(v), k) for k, v in by.items() if verdict(k) is None)
    if missing:
        for n, k in sorted(missing, reverse=True):
            print(f"UNMATCHED {n} {k}")
        if "--check" in sys.argv:
            sys.exit(1)
    clusters = []
    for sig, rs in sorted(by.items(), key=lambda kv: (kv[1][0]["port"], -len(kv[1]), kv[0])):
        v = verdict(sig)
        if v is None:
            continue
        i, (_, cls, evidence, rationale) = v
        ids = sorted(r["id"] for r in rs)
        clusters.append({
            "signature": sig, "port": rs[0]["port"], "group": rs[0]["group"], "class_outcome": rs[0]["class"],
            "detail": rs[0]["detail"], "feature": rs[0]["feature"], "n": len(rs),
            "ops": dict(collections.Counter(r["op"] for r in rs)),
            "sources": dict(collections.Counter(r["source"] for r in rs)),
            "class": cls, "rule": i, "evidence": evidence, "rationale": rationale,
            "examined": picks(ids), "case_ids": ids,
        })
    with (HERE / "clusters.jsonl").open("w") as fh:
        for c in clusters:
            fh.write(json.dumps(c) + "\n")

    _, expected, _, scored = load()
    summary = {}
    for port in PORTS:
        cnt = collections.Counter()
        held = collections.Counter()
        for c in clusters:
            if c["port"] == port:
                cnt[c["class"]] += c["n"]
                held[c["class"]] += c["sources"].get("heldout", 0)
        outs = collections.Counter(s["outcome"] for s in scored[port].values())
        outs_h = collections.Counter(s["outcome"] for s in scored[port].values() if s["source"] == "heldout")

        def fid(o: collections.Counter, k: collections.Counter) -> dict:
            supported = sum(o.values()) - o["unsupported"]
            denom = supported - sum(k[c] for c in EXCLUDED)
            return {"pass": o["pass"], "scored": sum(o.values()), "supported": supported,
                    "raw_fidelity": round(o["pass"] / sum(o.values()), 4),
                    "supported_fidelity": round(o["pass"] / supported, 4) if supported else None,
                    "adjusted_denominator": denom, "adjusted_fidelity": round(o["pass"] / denom, 4) if denom else None}

        summary[port] = {"classes": {c: cnt[c] for c in CLASSES}, "classes_heldout": {c: held[c] for c in CLASSES},
                         "outcomes": dict(outs), "all": fid(outs, cnt), "heldout": fid(outs_h, held),
                         "clusters": sum(1 for c in clusters if c["port"] == port)}
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    for port, s in summary.items():
        print(port, s["clusters"], s["classes"], "heldout adj", s["heldout"]["adjusted_fidelity"],
              "all adj", s["all"]["adjusted_fidelity"])


if __name__ == "__main__":
    main()
