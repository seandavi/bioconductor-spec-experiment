"""Schema-check cases, resolved and expected files that exist."""

from __future__ import annotations

import json
import sys

from jsonschema import Draft202012Validator

from .common import BUILD, CASES, POOL_FILES, SPEC, iter_jsonl, load_ops, ops_in_call


def validator(name: str) -> Draft202012Validator:
    return Draft202012Validator(json.loads((SPEC / name).read_text()))


def integral_floats(node, path="$"):
    """CASES.md forbids bare integral doubles in a call (they must be {"dbl": [..]})."""
    if isinstance(node, float) and node.is_integer():
        yield path
    elif isinstance(node, dict):
        for k, v in node.items():
            if k != "dbl":
                yield from integral_floats(v, f"{path}.{k}")
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from integral_floats(v, f"{path}[{i}]")


def main() -> int:
    bad = 0
    cv, pv = validator("case.schema.json"), validator("projection.schema.json")
    ops = load_ops()
    seen: set[str] = set()
    n = 0
    for stem in POOL_FILES:
        path = CASES / f"{stem}.jsonl"
        if not path.exists():
            continue
        for case in iter_jsonl(path):
            n += 1
            errs = [e.message[:160] for e in cv.iter_errors(case)]
            errs += [f"bare integral double at {p}" for p in integral_floats(case["call"])]
            errs += [f"op outside ops.json: {o}" for o in ops_in_call(case["call"]) if o not in ops]
            if case["id"] in seen:
                errs.append("duplicate id")
            seen.add(case["id"])
            if errs:
                bad += 1
                if bad <= 20:
                    print(f"FAIL case {case['id']}: {errs[:3]}")
    print(f"cases: {n} checked, {bad} invalid")
    exp = BUILD / "expected.jsonl"
    if exp.exists():
        m = b = 0
        for rec in iter_jsonl(exp):
            m += 1
            proj = rec.get("result")
            if proj is not None and (errs := [e.message[:160] for e in pv.iter_errors(proj)]):
                b += 1
                if b <= 20:
                    print(f"FAIL expected {rec['id']}: {errs[:2]}")
        print(f"expected: {m} checked, {b} invalid")
        bad += b
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
