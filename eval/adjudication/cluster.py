"""Cluster every adjudicated (mismatch / port_error_unexpected) outcome by a reproducible signature.

Signature = port | op family (or `*` for cross-op argument validation) | outcome class | decisive difference | decisive input feature.

- R errored, port returned a value: the R error message is bucketed (R_ERROR_BUCKETS).
- Port raised: the port's message is bucketed (PORT_ERROR_BUCKETS), else `<phase>:<type>`.
- R's value is itself defective (an R-side probe says so): feature `r_nonlocal_na`, whatever the
  port did.
- Values differ: `diagnose.decisive` names the representation-level normalization under which
  R and port agree (names, num_type, seqlevel_order, ...); when none does, the first diff kind/path
  plus the highest-priority input feature (INPUT_FEATURES) names the cluster.

Usage: python adjudication/cluster.py [--write]   (prints a summary; --write emits signatures.jsonl)
"""
from __future__ import annotations

import collections
import copy
import json
import re
import sys
from pathlib import Path

from diagnose import decisive, row_order
from load import ADJUDICATED, EVAL, PORTS, load

FAMILY = {o["op"]: o["family"] for o in json.loads((EVAL / "spec" / "ops.json").read_text())["ops"]}
# argument-validation buckets are one behaviour across every op: clustered with op "*"
CROSS_OP = {"flag_validation", "scalar_validation", "type_validation", "length_validation",
            "rejects_recycling", "rejects_double_arg", "na_arg_crash"}

HERE = Path(__file__).resolve().parent
OUT = HERE / "signatures.jsonl"
# R-side probe results (probe_*.py): case ids whose R value is itself defective
R_NONLOCAL_NA = {r["id"] for r in map(json.loads, (HERE / "probe_coverage_na.jsonl").open()) if r["local_na"] is False}
INT32_EDGE = 2**30

R_ERROR_BUCKETS: list[tuple[str, str]] = [
    ("flag_validation", r"must be TRUE or FALSE|must be TRUE of FALSE|must be a single logical"),
    ("scalar_validation", r"must be a single (integer|number|string|numeric)|must be non-negative|"
                          r"must be a non-negative|must be >= 0|cannot be NA|must be a single"),
    ("seqinfo_incompatible", r"incompatible seqlengths|not the same as the|incompatible"),
    ("int32_overflow", r"2\^31|integer overflow|invalid class .IRanges. object|NAs produced by integer overflow"),
    ("na_in_logic", r"missing value where TRUE/FALSE needed|logical subscript contains NAs"),
    ("na_input_rejected", r"cannot contain NAs|contains NAs|missing values not allowed"),
    ("rank_default_ties", r"'arg' must be of length 1"),
    ("length_validation", r"must have the same length|must have length|length of|recycl"),
    ("type_validation", r"must be (an|a) (integer|numeric|character)|must be a numeric|must be a character|"
                        r"cannot subset by character|bad supplied"),
]

PORT_ERROR_BUCKETS: list[tuple[str, str]] = [
    ("rejects_recycling", r"has length \d+ but \d+ was expected|must be a single value or have the length|"
                          r"must have length \d+ or the length|has \d+ rows but object has length|"
                          r"length of the supplied .* must be equal|lengths? .* (differ|must match)|"
                          r"not a multiple|must have the same length|recycl"),
    ("rejects_double_arg", r"whole numbers|must be integer|integer-valued|not integral|fraction"),
    ("int32_overflow", r"2\^31|overflow|int32"),
    ("mcols_column_name_arg", r"could not convert string to float|mcols column|column .* not found"),
    ("na_arg_crash", r"cannot convert float NaN to integer|NaN|masked"),
]

INPUT_FEATURES = ("int32_edge", "circular", "coord_le0", "na_arg", "empty", "zero_width", "strand",
                  "multi_seqlevel", "names", "mcols")


def _bucket(msg: str, buckets: list[tuple[str, str]]) -> str | None:
    for name, pat in buckets:
        if re.search(pat, msg):
            return name
    return None


def _numbers(node) -> list[int | float]:
    out: list = []
    if isinstance(node, dict):
        for v in node.values():
            out += _numbers(v)
    elif isinstance(node, list):
        for v in node:
            out += _numbers(v)
    elif isinstance(node, (int, float)) and not isinstance(node, bool):
        out.append(node)
    return out


def input_features(case: dict, exp: dict) -> set[str]:
    f = set(exp.get("features", []))
    if case["id"] in R_NONLOCAL_NA:
        f.add("r_nonlocal_na")
    if any(abs(v) >= INT32_EDGE for v in _numbers([case["inputs"], case["call"]])):
        f.add("int32_edge")
    for p in case["inputs"].values():
        if isinstance(p, dict) and p.get("kind") in ("IRanges", "GRanges"):
            if any(s is not None and (s <= 0 or (w is not None and s + w - 1 <= 0))
                   for s, w in zip(p["start"], p["width"])):
                f.add("coord_le0")
    return f


def _msg(m: str) -> str:
    return re.sub(r"\d+", "N", (m or "").split("\n")[0])[:60]


def signature(port: str, case: dict, exp: dict, rec: dict, scored: dict) -> dict:
    op, outcome = scored["op"], scored["outcome"]
    feats = input_features(case, exp)
    feature = next((f for f in INPUT_FEATURES if f in feats), "plain")
    if outcome == "port_error_unexpected":
        err = rec.get("error", {})
        b = _bucket(err.get("message", ""), PORT_ERROR_BUCKETS)
        if rec.get("phase") == "build" and b is None:
            b = "build"
        detail = b or f"{err.get('type')}:{_msg(err.get('message', ''))}"
        cls = f"port_error:{rec.get('phase')}"
        if err.get("origin") == "timeout":
            detail = "timeout"
    elif scored["rule"] == "row4":
        cls = "r_error_port_value"
        detail = _bucket(exp["result"]["message"], R_ERROR_BUCKETS) or f"other:{_msg(exp['result']['message'])}"
    elif scored["rule"] == "R3":
        cls, detail = "mismatch", "R3:" + ",".join(scored["detail"]["kinds"])
    else:
        r, p = exp["result"], rec["result"]
        if scored["rule"] == "R2":
            r, p = row_order(copy.deepcopy(r)), row_order(copy.deepcopy(p))
        cls, detail = "mismatch", decisive(r, p)
    if not detail.startswith(("diff:", "R3:")) and cls == "mismatch":
        feature = "-"  # a representation normalization explains it; the input feature is not decisive
    if cls != "mismatch":
        feature = "-"
    if "r_nonlocal_na" in feats:  # R's own value is defective here, whatever the port did (probe_coverage_na.py)
        feature = "r_nonlocal_na"
    if detail.startswith("diff:") or (detail.startswith("R3:") and feature != "plain"):
        detail = "diff"  # no normalization explains it: the input feature names the cluster
    group = "*" if detail in CROSS_OP else FAMILY[op]
    sig = f"{port}|{group}|{cls}|{detail}|{feature}"
    return {"signature": sig, "port": port, "group": group, "op": op, "domain": exp["domain"], "class": cls,
            "detail": detail, "feature": feature, "id": case["id"], "source": case["source"],
            "features": sorted(feats)}


def main() -> None:
    resolved, expected, raw, scored = load()
    rows = []
    for port in PORTS:
        for cid, s in scored[port].items():
            if s["outcome"] in ADJUDICATED:
                rows.append(signature(port, resolved[cid], expected[cid], raw[port][cid], s))
    if "--write" in sys.argv:
        with OUT.open("w") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
    by_port = collections.Counter(r["port"] for r in rows)
    sigs = collections.Counter((r["port"], r["signature"]) for r in rows)
    for port in PORTS:
        n = sum(1 for (p, _) in sigs if p == port)
        print(f"{port}: {by_port[port]} outcomes, {n} signatures")
    if "--list" in sys.argv:
        port = sys.argv[sys.argv.index("--list") + 1]
        for (p, s), n in sigs.most_common():
            if p == port:
                print(n, s)


if __name__ == "__main__":
    main()
