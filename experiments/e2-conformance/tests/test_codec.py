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


# ---- review finding 3: IntegerList / Factor columns are typed by content, not by column name -------


def _dframe(columns: dict[str, object], nrow: int) -> dict[str, object]:
    return {"$type": "DFrame", "nrow": nrow, "columns": columns}


@pytest.mark.parametrize("name", ["groups", "revmap"])
def test_integerlist_column_roundtrips_under_any_name(name: str) -> None:
    il = {"$type": "IntegerList", "elements": [[1, 2], [], [3]], "names": None}
    tree = _dframe({name: il}, 3)
    assert compare(encode(decode(tree)), tree)


def test_only_revmap_is_index_shifted() -> None:
    from ranges_e2.conformance.codec import DFrameOut

    il = {"$type": "IntegerList", "elements": [[1, 2]], "names": None}
    g = decode(_dframe({"groups": il, "revmap": il}, 1))
    assert isinstance(g, DFrameOut) and g.frame is not None
    assert g.frame["groups"].to_list() == [[1, 2]]  # plain values, as written
    assert g.frame["revmap"].to_list() == [[0, 1]]  # 0-based indices in the Python API


def test_factor_column_roundtrips_with_levels_and_na() -> None:
    f = {"$type": "Factor", "levels": {"$chr": ["lo", "hi", "unused"]}, "codes": {"$int": [2, None, 1]}}
    tree = _dframe({"f": f}, 3)
    assert compare(encode(decode(tree)), tree)
