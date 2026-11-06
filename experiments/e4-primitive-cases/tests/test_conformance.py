"""One pytest param per case: build natively, run, project, compare to R's golden."""

from __future__ import annotations

from typing import Any

import pytest

from conformance.build import build_inputs
from conformance.compare import compare
from conformance.loader import Case, load_cases
from conformance.ops import run_op
from conformance.project import project

CASES = load_cases()


def _run(spec: dict[str, Any]) -> dict[str, Any]:
    inputs = build_inputs(spec["inputs"])
    result, is_index = run_op(spec["call"], inputs)
    return project(result, is_index)


@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_case(case: Case) -> None:
    if case.golden is None:
        pytest.fail("no golden for this case: run `just goldens`")
    py = case.spec.get("python") or {}
    if "skip" in py:
        pytest.skip(py["skip"])
    try:
        actual = _run(case.spec)
    except (ValueError, TypeError) as e:
        if case.golden["kind"] == "error" or py.get("expect_error"):
            return
        raise AssertionError(f"Python raised {type(e).__name__}: {e}; R succeeded") from e
    if case.golden["kind"] == "error":
        pytest.fail(f"R raised ({case.golden['message']!r}); Python returned a result")
    if py.get("expect_error"):
        pytest.fail("python.expect_error set but Python succeeded")
    problems = compare(case.spec, case.golden, actual)
    assert not problems, "\n".join(problems)
