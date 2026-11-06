"""R probe: what does R return for each circular-feature case when every sequence is linear?

Every `isCircular` TRUE in the case's input literals is set to FALSE and the case is re-run in R.
`circular_matters` = R's result changes. A port whose value equals `r_linear` (and not R's) is
computing linear semantics on a circular sequence. Writes probe_circular.jsonl:
{id, circular_matters, r_linear}.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

from load import EVAL
from rprobe import literal_cases, run

OUT = Path(__file__).resolve().parent / "probe_circular.jsonl"


def linearize(node):
    if isinstance(node, dict):
        out = {}
        for k, v in node.items():
            if k == "isCircular" and isinstance(v, dict) and "lgl" in v:
                v = {**v, "lgl": [False if b else b for b in v["lgl"]]}
            out[k] = linearize(v)
        return out
    if isinstance(node, list):
        return [linearize(v) for v in node]
    return node


def main() -> None:
    expected = {r["id"]: r for r in map(json.loads, (EVAL / "build" / "expected.jsonl").open())}
    todo = []
    for cid, c in literal_cases().items():
        if "circular" in expected.get(cid, {}).get("features", []):
            lc = copy.deepcopy(c)
            lc["inputs"] = linearize(lc["inputs"])
            lc["id"] = cid + "#linear"
            todo.append(lc)
    got = run(todo)
    n = 0
    with OUT.open("w") as fh:
        for lc in todo:
            cid = lc["id"].removesuffix("#linear")
            lin = got[lc["id"]].get("result")
            matters = lin != expected[cid]["result"]
            n += matters
            fh.write(json.dumps({"id": cid, "circular_matters": matters, "r_linear": lin}) + "\n")
    print(f"{len(todo)} circular-feature cases; R's result depends on circularity in {n}")


if __name__ == "__main__":
    main()
