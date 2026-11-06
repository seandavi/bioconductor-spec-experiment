"""Check the contract artifacts against their examples.

Usage: uv run --no-project --with jsonschema python eval/spec/examples/validate.py
"""

import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator

SPEC = Path(__file__).resolve().parent.parent


def lines(name: str) -> list[dict]:
    text = (SPEC / "examples" / name).read_text()
    return [json.loads(line) for line in text.splitlines() if line.strip()]


def integral_floats(node, path="$"):
    """JSON Schema sees 2.0 as an integer; CASES.md forbids bare integral doubles."""
    if isinstance(node, float) and node.is_integer():
        yield path
    elif isinstance(node, dict):
        # typed vector payloads may hold integral doubles: {"dbl": [2.0]}
        for k, v in node.items():
            if k != "dbl":
                yield from integral_floats(v, f"{path}.{k}")
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from integral_floats(v, f"{path}[{i}]")


def case_errors(v: Draft202012Validator, case: dict) -> list[str]:
    errs = [e.message[:200] for e in v.iter_errors(case)]
    errs += [f"bare integral double at {p}" for p in integral_floats(case.get("call", {}))]
    return errs


def main() -> int:
    ok = True
    checks = [
        ("case.schema.json", "cases.jsonl", "invalid-cases.jsonl", case_errors),
        ("projection.schema.json", "projections.jsonl", "invalid.jsonl",
         lambda v, d: [e.message[:200] for e in v.iter_errors(d)]),
    ]
    for schema_name, good, bad, errors_of in checks:
        schema = json.loads((SPEC / schema_name).read_text())
        Draft202012Validator.check_schema(schema)
        v = Draft202012Validator(schema)
        for i, d in enumerate(lines(good)):
            if errs := errors_of(v, d):
                ok = False
                print(f"FAIL: {good}:{i + 1} rejected: {errs[:3]}")
        for i, d in enumerate(lines(bad)):
            if not errors_of(v, d):
                ok = False
                print(f"FAIL: {bad}:{i + 1} accepted")
        print(f"{schema_name}: {len(lines(good))} valid, {len(lines(bad))} invalid examples")

    ops = json.loads((SPEC / "ops.json").read_text())
    names = [o["op"] for o in ops["ops"]]
    families = set(ops["families"])
    assert len(names) == len(set(names)), "duplicate op"
    assert {o["family"] for o in ops["ops"]} == families, "family set mismatch"
    assert all(set(o["positional"]) <= set(o["named"]) for o in ops["ops"]), "positional not in named"

    def ops_in(call):
        yield call["op"]
        for a in [*call["args"], *call["kwargs"].values()]:
            if isinstance(a, dict) and "$call" in a:
                yield from ops_in(a["$call"])

    for d in lines("cases.jsonl"):
        unknown = [op for op in ops_in(d["call"]) if op not in names]
        if unknown:
            ok = False
            print(f"FAIL: {d['id']} uses ops outside ops.json: {unknown}")
    print(f"ops.json: {len(names)} ops, {len(families)} families")
    print("OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
