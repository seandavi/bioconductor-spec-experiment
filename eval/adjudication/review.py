"""Print the examined cases of every cluster no verdict rule matches yet (adjudication work list).

Usage: python adjudication/review.py PORT [SKIP_SUBSTRING ...]
"""
from __future__ import annotations

import collections
import json
import sys

from build import verdict
from cluster import OUT
from load import load
from show import _compact, _short, picks


def main() -> None:
    port, skip = sys.argv[1], sys.argv[2:]
    resolved, expected, raw, scored = load()
    by = collections.defaultdict(list)
    for r in map(json.loads, OUT.open()):
        if r["port"] == port and verdict(r["signature"]) is None and not any(s in r["signature"] for s in skip):
            by[r["signature"]].append(r["id"])
    for sig, ids in sorted(by.items(), key=lambda x: (-len(x[1]), x[0])):
        print(f"\n######## {len(ids)} {sig}")
        for cid in picks(ids):
            p = raw[port][cid]
            print("--", cid)
            print("  IN", _short(_compact(resolved[cid]), 500))
            print("  R ", _short(expected[cid]["result"], 300))
            print("  P ", _short(p.get("result") or p.get("error"), 300))


if __name__ == "__main__":
    main()
