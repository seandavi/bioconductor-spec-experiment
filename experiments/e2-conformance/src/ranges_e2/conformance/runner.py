"""Conformance runner: load cases + R goldens, run every case against ranges_e2, report.

    uv run python -m ranges_e2.conformance.runner [--report PATH] [group ...]
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .codec import compare, decode, encode
from .ops import run_op

CONFORMANCE_DIR = Path(__file__).resolve().parents[3] / "conformance"
ALLOWED_ERRORS = (ValueError, TypeError, IndexError, KeyError)


@dataclass
class Case:
    id: str
    group: str
    clause: str
    op: str
    args: list[Any]
    kwargs: dict[str, Any]
    compare: str
    fixtures: dict[str, Any] = field(repr=False)


def load_cases(base: Path = CONFORMANCE_DIR, groups: list[str] | None = None) -> list[Case]:
    cases: list[Case] = []
    for path in sorted((base / "cases").glob("*.json")):
        doc = json.loads(path.read_text())
        if groups and doc["group"] not in groups:
            continue
        for c in doc["cases"]:
            cases.append(
                Case(
                    c["id"],
                    doc["group"],
                    c["clause"],
                    c["op"],
                    c.get("args", []),
                    c.get("kwargs", {}),
                    c.get("compare", doc.get("compare", "exact")),
                    doc["fixtures"],
                )
            )
    return cases


def validate(base: Path = CONFORMANCE_DIR) -> list[str]:
    """Same rules as the R runner: unique ids, declared clauses, every clause covered, refs resolve."""
    errs: list[str] = []
    seen: set[str] = set()
    for path in sorted((base / "cases").glob("*.json")):
        doc = json.loads(path.read_text())
        used: set[str] = set()
        for c in doc["cases"]:
            if c["id"] in seen:
                errs.append(f"duplicate id {c['id']}")
            seen.add(c["id"])
            if c["clause"] not in doc["clauses"]:
                errs.append(f"{c['id']}: unknown clause {c['clause']}")
            used.add(c["clause"])
            for ref in _refs([c.get("args"), c.get("kwargs")]):
                if ref not in doc["fixtures"]:
                    errs.append(f"{c['id']}: unresolved ref {ref}")
        errs += [f"{doc['group']}: clause without case {cl}" for cl in doc["clauses"] if cl not in used]
    return errs


def _refs(node: Any) -> list[str]:
    if isinstance(node, dict):
        own = [node["$ref"]] if "$ref" in node else []
        return own + [r for v in node.values() for r in _refs(v)]
    if isinstance(node, list):
        return [r for v in node for r in _refs(v)]
    return []


def load_expected(base: Path = CONFORMANCE_DIR) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for path in sorted((base / "expected").glob("*.json")):
        out.update(json.loads(path.read_text())["cases"])
    return out


def load_xfail(base: Path = CONFORMANCE_DIR) -> dict[str, str]:
    path = base / "xfail.json"
    return json.loads(path.read_text()) if path.exists() else {}


_fixture_cache: dict[str, dict[str, Any]] = {}


def decode_fixtures(case: Case) -> dict[str, Any]:
    key = case.group
    if key not in _fixture_cache:
        fx: dict[str, Any] = {}
        for name, tree in case.fixtures.items():
            fx[name] = decode(tree, fx)
        _fixture_cache[key] = fx
    return _fixture_cache[key]


def run_case(case: Case, expected: dict[str, Any]) -> tuple[bool, str]:
    """Return (matches R golden, detail). Raises on harness errors (decode/encode failures)."""
    entry = expected[case.id]
    fx = decode_fixtures(case)
    args = [decode(a, fx) for a in case.args]
    kwargs = {k: decode(v, fx) for k, v in case.kwargs.items()}
    if "error" in entry:
        try:
            run_op(case.op, args, kwargs)
        except ALLOWED_ERRORS:
            return True, ""
        except NotImplementedError:
            raise
        return False, "R raised an error; ranges_e2 returned a value"
    try:
        result = run_op(case.op, args, kwargs)
    except ALLOWED_ERRORS as e:
        return False, f"ranges_e2 raised {type(e).__name__}: {e}"
    ok = compare(encode(result), entry["value"], case.compare)
    return ok, "" if ok else "value differs from R golden"


def status_of(case: Case, expected: dict[str, Any], xfail: dict[str, str]) -> tuple[str, str]:
    try:
        ok, detail = run_case(case, expected)
        status = "pass" if ok else "fail"
    except Exception as e:  # noqa: BLE001 - any crash is a case failure, reported with its type
        status, detail = "error", f"{type(e).__name__}: {e}"
    if case.id in xfail:
        status = {"pass": "xpass", "fail": "xfail", "error": "xfail"}[status]
    return status, detail


def run_all(groups: list[str] | None = None) -> dict[str, tuple[str, str]]:
    expected = load_expected()
    xfail = load_xfail()
    return {c.id: status_of(c, expected, xfail) for c in load_cases(groups=groups)}


def main(argv: list[str]) -> int:
    report_path = None
    groups: list[str] = []
    it = iter(argv)
    for a in it:
        if a == "--report":
            report_path = Path(next(it))
        else:
            groups.append(a)
    errs = validate()
    if errs:
        print("INVALID:", *errs, sep="\n  ")
        return 2
    results = run_all(groups or None)
    counts: dict[str, int] = {}
    for cid, (st, detail) in results.items():
        counts[st] = counts.get(st, 0) + 1
        if st in ("fail", "error", "xpass"):
            print(f"{st.upper():5} {cid} {detail}")
    print("python:", ", ".join(f"{k}={v}" for k, v in sorted(counts.items())), f"total={len(results)}")
    if report_path and not groups:
        report_path.write_text(
            json.dumps({"implementation": "python", "cases": {k: v[0] for k, v in results.items()}}, indent=1) + "\n"
        )
    return 0 if not {"fail", "error", "xpass"} & counts.keys() else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
