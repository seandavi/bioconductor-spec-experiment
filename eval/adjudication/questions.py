"""Q4 tables: circular semantics per port and family (from probe_circular.jsonl), written to q4.json.

For every case where R's result depends on circularity (circular_matters), each port's outcome is
`matches_R` (pass), `matches_R_linear` (equals R's answer with all sequences linear), `unsupported`,
or `other`. Accessor/constructor/seqinfo/metadata/subset families are left out (circularity there is
just a stored flag).
"""
from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

from load import EVAL, PORTS, load

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(EVAL / "src"))
from cluster import FAMILY, _equals_linear  # noqa: E402

SKIP = {"seqinfo", "accessors", "construct", "metadata", "subset_concat"}


def main() -> None:
    _, expected, raw, scored = load()
    probe = [r for r in map(json.loads, (HERE / "probe_circular.jsonl").open()) if r["circular_matters"]]
    table: dict[str, dict[str, collections.Counter]] = {p: collections.defaultdict(collections.Counter) for p in PORTS}
    for r in probe:
        fam = FAMILY[expected[r["id"]]["op"]]
        if fam in SKIP:
            continue
        for p in PORTS:
            o = scored[p][r["id"]]["outcome"]
            rec = raw[p][r["id"]]
            if o == "pass":
                k = "matches_R"
            elif o == "unsupported":
                k = "unsupported"
            elif rec.get("status") == "value" and _equals_linear(r["id"], rec["result"]):
                k = "matches_R_linear"
            else:
                k = "other"
            table[p][fam][k] += 1
    out = {"circular": {p: {f: dict(c) for f, c in sorted(fams.items())} for p, fams in table.items()},
           "seqlevels": json.loads((HERE / "probe_seqlevels.json").read_text())}
    (HERE / "q4.json").write_text(json.dumps(out, indent=1) + "\n")
    for p in PORTS:
        print(p, {f: dict(c) for f, c in sorted(table[p].items())})


if __name__ == "__main__":
    main()
