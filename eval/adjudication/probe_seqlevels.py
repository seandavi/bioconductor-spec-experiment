"""Seqlevel order of a GRanges built without explicit seqinfo: R vs every port (Q4).

Runs a few literal cases through R (oracle/oracle.R) and through each port's adapter
(runner/run_port.py in the port's own venv; ports are only read). Writes probe_seqlevels.json:
{case_id: {"call": ..., "R": seqlevels, "<port>": seqlevels | error}}.
"""
from __future__ import annotations

import json
import subprocess
import tempfile
import tomllib
from pathlib import Path

from load import EVAL, PORTS
from rprobe import run

OUT = Path(__file__).resolve().parent / "probe_seqlevels.json"


def _ir(n: int) -> dict:
    return {"IRanges": {"start": {"int": list(range(1, n + 1))}, "width": {"int": [1] * n}}}


CASES = [
    {"id": "chr-mixed", "inputs": {"r": _ir(6)},
     "call": {"op": "GRanges", "args": [], "kwargs": {
         "seqnames": {"chr": ["chr2", "chr10", "chr1", "chr2", "chrX", "chr1"]}, "ranges": {"$ref": "r"}}}},
    {"id": "letters", "inputs": {"r": _ir(4)},
     "call": {"op": "GRanges", "args": [], "kwargs": {
         "seqnames": {"chr": ["b", "a", "c", "a"]}, "ranges": {"$ref": "r"}}}},
    {"id": "factor-levels", "inputs": {"r": _ir(3), "f": {"Factor": {"levels": ["z", "a", "q"], "codes": [2, 1, 2]}}},
     "call": {"op": "GRanges", "args": [], "kwargs": {"seqnames": {"$ref": "f"}, "ranges": {"$ref": "r"}}}},
    {"id": "seqlengths-names", "inputs": {"r": _ir(2)},
     "call": {"op": "GRanges", "args": [], "kwargs": {
         "seqnames": {"chr": ["a", "b"]}, "ranges": {"$ref": "r"},
         "seqlengths": {"int": [100, 200], "names": ["b", "a"]}}}},
    {"id": "c-merge", "inputs": {
        "x": {"GRanges": {"seqnames": {"chr": ["chr2"]}, "start": {"int": [1]}, "width": {"int": [2]}}},
        "y": {"GRanges": {"seqnames": {"chr": ["chr1", "chr10"]}, "start": {"int": [1, 3]}, "width": {"int": [2, 2]}}}},
     "call": {"op": "c", "args": [{"$ref": "x"}, {"$ref": "y"}], "kwargs": {}}},
]


def main() -> None:
    r = run(CASES)
    resolved = [{"id": c["id"], "source": "probe", "op": c["call"]["op"], "inputs": r[c["id"]]["inputs"],
                 "call": c["call"]} for c in CASES]
    table = {c["id"]: {"call": c["call"], "R": r[c["id"]]["result"].get("seqinfo", {}).get("seqnames")
                       if r[c["id"]]["result"]["kind"] == "GRanges" else r[c["id"]]["result"]} for c in CASES}
    cfg = tomllib.loads((EVAL / "ports.toml").read_text())
    with tempfile.TemporaryDirectory() as d:
        src = Path(d) / "resolved.jsonl"
        src.write_text("".join(json.dumps(x) + "\n" for x in resolved))
        for port in PORTS:
            p = cfg["ports"][port]["project"]
            proj = EVAL.parent / p if p.startswith("eval/") else Path(cfg["WORKTREES"]) / p
            out = Path(d) / f"{port}.jsonl"
            subprocess.run(["uv", "run", "--frozen", "--project", str(proj), "python",
                            str(EVAL / "runner" / "run_port.py"), "--port", port, "--resolved", str(src),
                            "--out", str(out)], check=True, cwd=EVAL.parent,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for rec in map(json.loads, out.open()):
                res = rec.get("result") or {}
                table[rec["id"]][port] = (res.get("seqinfo", {}).get("seqnames") if res.get("kind") == "GRanges"
                                          else rec.get("error") or rec.get("reason") or rec.get("status"))
    OUT.write_text(json.dumps(table, indent=1) + "\n")
    for cid, row in table.items():
        print(cid, {k: v for k, v in row.items() if k != "call"})


if __name__ == "__main__":
    main()
