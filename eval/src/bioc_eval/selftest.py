"""Self-test verdicts: r_echo must score 100 % pass; mutant's outcome must equal its prediction on every case."""

from __future__ import annotations

import sys
from collections import Counter

from .common import BUILD, RESULTS, iter_jsonl


def main() -> int:
    ok = True
    echo = list(iter_jsonl(RESULTS / "scored" / "r_echo.jsonl"))
    c = Counter(r["outcome"] for r in echo)
    echo_ok = c["pass"] == len(echo) and len(echo) > 0
    print(f"r_echo: {c['pass']}/{len(echo)} pass {'OK' if echo_ok else 'FAIL ' + str(dict(c))}")
    ok &= echo_ok
    pred = {r["id"]: r["predicted"] for r in iter_jsonl(BUILD / "results" / "mutant.predicted.jsonl")}
    scored = list(iter_jsonl(RESULTS / "scored" / "mutant.jsonl"))
    wrong = [(r["id"], pred.get(r["id"]), r["outcome"], r.get("rule")) for r in scored if pred.get(r["id"]) != r["outcome"]]
    by = Counter(r["outcome"] for r in scored)
    perturbed = sum(1 for r in scored if pred.get(r["id"]) != "pass") + 0
    print(f"mutant: {len(scored)} cases, outcomes {dict(by)}, predicted-vs-scored disagreements: {len(wrong)}")
    for w in wrong[:10]:
        print("   ", w)
    ok &= not wrong and len(scored) > 0 and len(pred) == len(scored)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
