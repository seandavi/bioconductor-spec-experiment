"""decode(golden) -> encode must reproduce every R golden: the Python codec agrees with R's."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ranges_e2.conformance.codec import compare, decode, encode

EXPECTED = Path(__file__).resolve().parents[1] / "conformance" / "expected"


def _goldens() -> list[tuple[str, object]]:
    out = []
    for path in sorted(EXPECTED.glob("*.json")):
        for case_id, entry in json.loads(path.read_text())["cases"].items():
            if "value" in entry:
                out.append((case_id, entry["value"]))
    return out


GOLDENS = _goldens()


@pytest.mark.parametrize("case_id,value", GOLDENS, ids=[g[0] for g in GOLDENS])
def test_golden_roundtrips(case_id: str, value: object) -> None:
    assert compare(encode(decode(value)), value), case_id


def test_goldens_exist() -> None:
    assert len(GOLDENS) > 500
