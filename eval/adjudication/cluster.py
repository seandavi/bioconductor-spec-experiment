"""Cluster every adjudicated (mismatch / port_error_unexpected) outcome by a reproducible signature.

Signature = port | op family (or `*` for cross-op argument validation) | outcome class | decisive difference | decisive input feature.

- R errored, port returned a value: the R error message is bucketed (R_ERROR_BUCKETS).
- Port raised: the port's message is bucketed (PORT_ERROR_BUCKETS), else `<phase>:<type>`.
- R's value is itself defective (an R-side probe says so): feature `r_nonlocal_na`, whatever the
  port did.
- The port's value equals R's answer for the same case with all sequences linear
  (probe_circular.py): detail `equals_R_linear`.
- R's value changes when every range is shifted by +1000 (probe_translate.py): feature
  `r_translation_dependent`; coord_le0 cases R answers consistently become `coord_le0_r_invariant`.
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

from diagnose import decisive, diffs, row_order
from load import ADJUDICATED, EVAL, PORTS, load

FAMILY = {o["op"]: o["family"] for o in json.loads((EVAL / "spec" / "ops.json").read_text())["ops"]}
# argument-validation buckets are one behaviour across every op: clustered with op "*"
CROSS_OP = {"flag_validation", "scalar_validation", "type_validation", "length_validation",
            "rejects_recycling", "rejects_names_padding", "rejects_double_arg", "na_arg_crash"}

HERE = Path(__file__).resolve().parent
OUT = HERE / "signatures.jsonl"
# R-side probe results (probe_*.py): case ids whose R value is itself defective
R_NONLOCAL_NA = {r["id"] for r in map(json.loads, (HERE / "probe_coverage_na.jsonl").open()) if r["local_na"] is False}
# R's result for the same case with every sequence made linear (probe_circular.py)
# R result changes when every range is translated by +1000 (probe_translate.py): True / False / None
R_TRANSLATION = {r["id"]: r["translation_dependent"] for r in map(json.loads, (HERE / "probe_translate.jsonl").open())}
R_LINEAR = {r["id"]: r["r_linear"] for r in map(json.loads, (HERE / "probe_circular.jsonl").open())
            if r["circular_matters"]}
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
    ("rejects_names_padding", r"'names' has length \d+ but \d+ was expected"),
    ("rejects_scalar_genome", r"genome has length 1, expected"),
    ("rejects_names_padding", r"names has length \d+, expected \d+"),
    # pintersect(resolve.empty=<invalid>) where no empty intersection arises; 'start' is a valid partial match
    ("resolve_empty_validation", r"unknown resolve_empty (?!'start')|'?resolve_empty'? must be"),
    ("mcols_type_conflict_crash", r"SchemaError|cannot combine: mcols columns differ"),
    ("rejects_recycling", r"has length \d+ but \d+ was expected|must be a single value or have the length|"
                          r"must have length \d+ or \d+, got \d+|incompatible lengths: \w+=|requires equal lengths|"
                          r"has \d+ rows, expected \d+|must be scalar or one per range|"
                          r"has length \d+, expected \d+( or \d+)?|equal lengths or one of length|lengths do not match|"
                          r"need equal lengths|"
                          r"must have length \d+ or the length|has \d+ rows but object has length|"
                          r"lengths? (.* )?(differ|must match)|"
                          r"not a multiple|must have the same length|recycl"),
    ("rejects_double_arg", r"whole numbers|must be integer|integer-valued|not integral|fraction|non-integer value"),
    ("int32_overflow", r"2\^31|overflow|int32|integer bounds"),
    ("mcols_column_name_arg", r"could not convert string to float|mcols column|column .* not found"),
    ("na_arg_crash", r"cannot convert float NaN to integer|NaN|masked"),
]

INPUT_FEATURES = ("int32_edge", "circular", "coord_le0", "na_arg", "empty", "zero_width", "strand",
                  "multi_seqlevel", "names", "mcols")
# nearest-family R defects come from coordinates <= 0 (sentinel number line), so that feature leads there
NEAREST_FEATURES = ("int32_edge", "coord_le0", "empty", "circular") + INPUT_FEATURES[3:]
# GRanges set operations clip to [1, seqlength]: a range outside those bounds is the decisive input
SETOPS_FEATURES = ("out_of_bounds",) + INPUT_FEATURES


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
        if isinstance(p, dict) and p.get("kind") == "GRanges":
            lens = dict(zip(p["seqinfo"]["seqnames"], p["seqinfo"]["seqlengths"]))
            if any(s is not None and (s < 1 or (lens.get(q) is not None and s + w - 1 > lens[q]))
                   for q, s, w in zip(p["seqnames"], p["start"], p["width"])):
                f.add("out_of_bounds")
        if isinstance(p, dict) and p.get("kind") in ("IRanges", "GRanges"):
            if any(s is not None and (s <= 0 or (w is not None and s + w - 1 <= 0))
                   for s, w in zip(p["start"], p["width"])):
                f.add("coord_le0")
    return f


def _no_circ(o):
    """Projection with every Seqinfo isCircular dropped (the linearized R input changed it)."""
    o = row_order(copy.deepcopy(o))

    def walk(n):
        if isinstance(n, dict):
            if n.get("kind") == "Seqinfo":
                n["isCircular"] = None
            for v in n.values():
                walk(v)
        elif isinstance(n, list):
            for v in n:
                walk(v)
    walk(o)
    return o


def _equals_linear(cid: str, p) -> bool:
    lin = R_LINEAR.get(cid)
    return lin is not None and lin.get("kind") != "error" and not diffs(_no_circ(lin), _no_circ(p))


def _valid_resolve_empty(v) -> bool:
    """A resolve.empty value R's match.arg accepts: one of the choices or a unique prefix of one."""
    if isinstance(v, dict):
        v = v.get("chr", [None])[0] if len(v.get("chr", [])) == 1 else None
    choices = ("none", "max.start", "start.x")
    return isinstance(v, str) and sum(c.startswith(v) for c in choices) == 1


def _msg(m: str) -> str:
    return re.sub(r"\d+", "N", (m or "").split("\n")[0])[:60]


def signature(port: str, case: dict, exp: dict, rec: dict, scored: dict) -> dict:
    op, outcome = scored["op"], scored["outcome"]
    feats = input_features(case, exp)
    order = {"nearest": NEAREST_FEATURES, "setops": SETOPS_FEATURES}.get(FAMILY[op], INPUT_FEATURES)
    feature = next((f for f in order if f in feats), "plain")
    if outcome == "port_error_unexpected":
        err = rec.get("error", {})
        b = _bucket(f"{err.get('type')}: {err.get('message', '')}", PORT_ERROR_BUCKETS)
        if rec.get("phase") == "build" and b is None:
            b = "build"
        detail = b or f"{err.get('type')}:{_msg(err.get('message', ''))}"
        cls = f"port_error:{rec.get('phase')}"
        if detail == "resolve_empty_validation" and _valid_resolve_empty(case["call"]["kwargs"].get("resolve.empty")):
            detail = "rejects_valid_resolve_empty"  # R accepts this value (or its unique partial match) anyway
        if err.get("origin") == "timeout":
            detail = "timeout"
    elif scored["rule"] == "row4":
        cls = "r_error_port_value"
        detail = _bucket(exp["result"]["message"], R_ERROR_BUCKETS) or f"other:{_msg(exp['result']['message'])}"
    elif scored["rule"] == "R3":
        cls, detail = "mismatch", "R3:" + ",".join(scored["detail"]["kinds"])
        if _equals_linear(case["id"], rec["result"]):
            detail = "equals_R_linear"
    else:
        r, p = exp["result"], rec["result"]
        if scored["rule"] == "R2":
            r, p = row_order(copy.deepcopy(r)), row_order(copy.deepcopy(p))
        cls, detail = "mismatch", decisive(r, p)
        if _equals_linear(case["id"], rec["result"]):
            detail = "equals_R_linear"  # the port computed R's answer for linear sequences

    if detail == "equals_R_linear":
        feature = "circular"
    elif not detail.startswith(("diff:", "R3:")) and cls == "mismatch":
        feature = "-"  # a representation normalization explains it; the input feature is not decisive
    if cls.startswith("port_error"):
        feature = "-"
    w = case["call"]["kwargs"].get("weight") if op == "coverage" else None
    if feature == "-" and isinstance(w, dict) and any(v is None for k in ("int", "dbl") for v in w.get(k, [])):
        feature = "na_weight_int" if "int" in w else "na_weight_dbl"  # NA weights: storage type decides
    if R_TRANSLATION.get(case["id"]) is True:  # R's value moves with the coordinates (probe_translate.py)
        feature = "r_translation_dependent"
    elif feature == "coord_le0" and R_TRANSLATION.get(case["id"]) is False:
        feature = "coord_le0_r_invariant"  # R gives the same answer shifted by +1000: not the sentinel defect
    if "r_nonlocal_na" in feats:  # R's own value is defective here, whatever the port did (probe_coverage_na.py)
        feature = "r_nonlocal_na"
    if detail.startswith("diff:") or (detail.startswith("R3:") and feature != "plain"):
        detail = "diff"  # no normalization explains it: the input feature names the cluster
        kw = case["call"]["kwargs"]
        if op == "reduce" and kw.get("min.gapwidth") == 0:
            feature = "min_gapwidth_0"  # R's empty-range rules differ by min.gapwidth / drop.empty.ranges
        elif op == "reduce" and kw.get("drop.empty.ranges") is True:
            feature = "drop_empty_ranges"
        elif op == "distance" and feature == "int32_edge":
            feature = "int32_edge_distance"  # R documents NA for distances beyond int32; no sentinel involved
    if cls == "r_error_port_value":
        # an NA flag R rejects vs an R failure on valid input: the input decides
        feature = "-" if detail not in ("na_in_logic", "length_validation") else next(
            (f for f in ("na_arg", "int32_edge", "circular", "coord_le0") if f in feats), "-")
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
