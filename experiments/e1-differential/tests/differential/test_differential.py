"""One pytest item per differential case: R (cached) vs Python."""

from __future__ import annotations

import json

import pytest

import harness

ALL = harness.load_cases()
VISIBLE = [c for c in ALL if "#" not in c.id]


@pytest.mark.parametrize("case", VISIBLE, ids=[c.id for c in VISIBLE])
def test_case(case: harness.Case, request: pytest.FixtureRequest) -> None:
    force = bool(request.config.getoption("--refresh-r"))
    r = harness.r_result(case, ALL, force=force)
    if r is None:
        pytest.skip("R result stale/missing and R unavailable")
    universe = None
    if case.spec.get("compare") == "arbitrary":
        universe = harness.r_result(next(c for c in ALL if c.id == case.id + "#universe"), ALL)
    py = harness.py_result(case)
    fail = harness.evaluate(case, r, py, universe)
    if fail:
        pytest.fail(
            f"{fail}\n  R : {json.dumps(r.get('value', r.get('error')))[:1500]}"
            f"\n  Py: {json.dumps(py.get('value', py.get('error')))[:1500]}",
            pytrace=False,
        )
