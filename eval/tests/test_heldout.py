"""Heldout generator: deterministic, schema-valid, covers every op family and feature, no foreign data."""

import json
import re
import subprocess
import sys

from bioc_eval.common import CASES, EVAL, load_ops, read_jsonl


def test_generator_is_deterministic_and_committed_output_current():
    r = subprocess.run([sys.executable, str(EVAL / "heldout" / "generate.py"), "--check"], capture_output=True, text=True,
                       env={"PYTHONHASHSEED": "7", "PATH": "/usr/bin:/bin"})
    assert r.returncode == 0, r.stdout + r.stderr


def test_every_op_in_both_classes_where_applicable():
    cases = read_jsonl(CASES / "heldout.jsonl")
    seen = {(c["op"], c["id"].split("/")[2]) for c in cases}
    ops = set(load_ops())
    assert {op for op, _ in seen} == ops
    for op in ("findOverlaps", "reduce", "shift", "c", "[", "order", "coverage"):
        assert (op, "IRanges") in seen and (op, "GRanges") in seen


def test_ids_source_and_seed_namespace():
    cases = read_jsonl(CASES / "heldout.jsonl")
    assert len({c["id"] for c in cases}) == len(cases)
    assert all(c["source"] == "heldout" and re.match(r"^heldout/[^/]+/[^/]+/\d+$", c["id"]) for c in cases)
    seeds = (EVAL / "heldout" / "seeds.txt").read_text().split()
    assert len(set(seeds)) == len(seeds) > 3000
    sys.path.insert(0, str(EVAL / "heldout"))
    from gen_core import seed_for
    assert seed_for("findOverlaps", 0) == int(seeds[[int(s) for s in seeds].index(seed_for("findOverlaps", 0))])


def test_required_features_present_in_every_family():
    cov = (EVAL / "heldout" / "COVERAGE.md").read_text()
    rows = [l for l in cov.splitlines() if l.startswith("| ") and not l.startswith("| family") and l.count("|") > 15]
    assert len(rows) == 12
    for row in rows:
        cells = [c.strip() for c in row.strip("|").split("|")]
        assert all(int(c) >= 5 for c in cells[2:9]), row  # empty zero_width names mcols strand seqinfo na_arg


def test_no_expected_values_leak_into_cases():
    for case in read_jsonl(CASES / "heldout.jsonl")[:50]:
        assert set(case) == {"id", "source", "origin", "op", "inputs", "call"}
