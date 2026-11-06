"""Show a case's R and Python results side by side.

    uv run python tests/differential/show.py <case-id> [--rerun-r]
"""

from __future__ import annotations

import json
import sys

import harness


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    cid, force = sys.argv[1], "--rerun-r" in sys.argv
    allc = harness.load_cases()
    case = next((c for c in allc if c.id == cid), None)
    if case is None:
        print(f"no such case: {cid}")
        return 2
    r = harness.r_result(case, allc, force=force)
    py = harness.py_result(case)
    print("R  :", json.dumps(r, indent=1) if r else "(stale; R unavailable)")
    print("Py :", json.dumps(py, indent=1))
    if r:
        u = harness.r_result(next((c for c in allc if c.id == cid + "#universe"), case), allc)
        print("verdict:", harness.evaluate(case, r, py, u) or "PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
