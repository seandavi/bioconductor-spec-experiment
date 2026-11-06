"""Run variant cases through the eval's own oracle (oracle/oracle.R) to test R for self-consistency.

`run(cases)` takes case dicts in literal form (cases/*.jsonl: id, inputs, call) and returns
{id: oracle record}. The variant ids carry a suffix naming the transformation.
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

from load import EVAL

ORACLE_R = EVAL / "oracle" / "oracle.R"


def literal_cases() -> dict[str, dict]:
    out: dict[str, dict] = {}
    for f in ("e1", "e2", "e4", "heldout"):
        with (EVAL / "cases" / f"{f}.jsonl").open() as fh:
            for line in fh:
                c = json.loads(line)
                out[c["id"]] = c
    return out


def run(cases: list[dict]) -> dict[str, dict]:
    with tempfile.TemporaryDirectory() as d:
        src, dst = Path(d) / "in.jsonl", Path(d) / "out.jsonl"
        with src.open("w") as fh:
            for c in cases:
                fh.write(json.dumps({"id": c["id"], "inputs": c["inputs"], "call": c["call"],
                                     "universe_call": None}) + "\n")
        subprocess.run(["Rscript", str(ORACLE_R), str(src), str(dst)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        with dst.open() as fh:
            return {r["id"]: r for r in map(json.loads, fh)}
