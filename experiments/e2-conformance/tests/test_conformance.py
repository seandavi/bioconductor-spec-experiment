"""Every conformance case, run against ranges_e2 and compared to the R golden.

`conformance/xfail.json` entries are strict xfails: they must keep failing (an unexpected pass
fails the run, forcing the entry and its NOTES record to be reviewed).
"""

from __future__ import annotations

import pytest

from ranges_e2.conformance import runner

CASES = runner.load_cases()
EXPECTED = runner.load_expected()
XFAIL = runner.load_xfail()


def _param(case: runner.Case) -> object:
    marks = []
    if case.id in XFAIL:
        marks.append(pytest.mark.xfail(reason=XFAIL[case.id], strict=True))
    return pytest.param(case, id=case.id, marks=marks)


def test_suite_is_valid() -> None:
    assert runner.validate() == []


def test_every_case_has_a_golden() -> None:
    assert {c.id for c in CASES} == set(EXPECTED)


def test_xfail_entries_name_real_cases() -> None:
    assert set(XFAIL) <= {c.id for c in CASES}


@pytest.mark.parametrize("case", [_param(c) for c in CASES])
def test_case(case: runner.Case) -> None:
    ok, detail = runner.run_case(case, EXPECTED)
    assert ok, f"{case.clause}: {detail}"
