"""Self-tests for the differential harness (no R needed)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "differential"))
import codec  # noqa: E402
import harness  # noqa: E402


def test_first_diff_numeric_leniency_and_strictness() -> None:
    a = {"class": "vector", "type": "int", "values": [1, 2]}
    b = {"class": "vector", "type": "num", "values": [1.0, 2.0]}
    assert codec.first_diff(a, b) is None
    assert codec.first_diff(a, {**b, "values": [1.0, 2.5]}) is not None
    assert codec.first_diff({"type": "lgl"}, {"type": "chr"}) is not None
    assert codec.first_diff({"error": "x"}, {"error": "completely different"}) is None
    assert codec.first_diff({"error": "x"}, {"v": 1}) is not None
    assert codec.first_diff([1, None], [1, None]) is None
    assert codec.first_diff([1], [1, 2]) is not None


def test_mcols_column_order_matters_other_dict_order_does_not() -> None:
    assert codec.first_diff({"mcols": {"a": 1, "b": 2}}, {"mcols": {"b": 2, "a": 1}}) is not None
    assert codec.first_diff({"x": 1, "y": 2}, {"y": 2, "x": 1}) is None


def test_case_ids_are_unique_and_files_parse() -> None:
    cases = harness.load_cases()
    ids = [c.id for c in cases]
    assert len(ids) == len(set(ids)) and len(ids) > 500


def test_case_hash_tracks_inputs_not_tags() -> None:
    base = harness.Case("c", "f", {"id": "c", "inputs": {"x": 1}, "call": {"fn": "f"}, "tags": ["a"]})
    retagged = harness.Case("c", "f", {"id": "c", "inputs": {"x": 1}, "call": {"fn": "f"}, "tags": ["b"]})
    changed = harness.Case("c", "f", {"id": "c", "inputs": {"x": 2}, "call": {"fn": "f"}})
    assert base.hash == retagged.hash != changed.hash


def test_every_cached_case_has_an_r_result() -> None:
    for c in harness.load_cases():
        doc = harness._load_cache(c.file)  # noqa: SLF001
        assert c.id in doc["entries"], f"{c.id}: run pytest with R available to fill the cache"


def test_decode_vectors() -> None:
    assert codec.decode([1, 2]).tolist() == [1, 2]
    assert codec.decode([1.5, 2]).dtype.kind == "f"
    assert codec.decode([True, False]).dtype == bool
    assert codec.decode([1, None]).mask.tolist() == [False, True]
    assert codec.decode({"chr": []}) == []
    assert codec.decode({"int": []}).dtype.kind == "i"
    assert codec.decode({"NULL": []}) is None


def test_arbitrary_comparison_accepts_any_member_of_the_hit_set() -> None:
    universe = {"value": {"class": "Hits", "from": [1, 1, 3], "to": [2, 5, 4], "nLnode": 3, "nRnode": 5}}
    ok = {"value": {"class": "vector", "type": "int", "values": [5, None, 4]}}
    bad_member = {"value": {"class": "vector", "type": "int", "values": [3, None, 4]}}
    missing = {"value": {"class": "vector", "type": "int", "values": [2, 1, 4]}}
    assert harness._arbitrary_ok(universe, ok) is None  # noqa: SLF001
    assert harness._arbitrary_ok(universe, bad_member) is not None  # noqa: SLF001
    assert harness._arbitrary_ok(universe, missing) is not None  # noqa: SLF001


def test_divergence_pins_python_value() -> None:
    case = harness.Case("d", "f", {"id": "d", "divergence": {"adr": "0000", "python": {"error": True}}})
    assert harness.evaluate(case, {"value": 1}, {"error": "x"}) is None
    assert harness.evaluate(case, {"value": 1}, {"value": 1}) is not None


def test_python_internal_errors_fail_even_if_r_errors() -> None:
    case = harness.Case("e", "f", {"id": "e"})
    msg = harness.evaluate(case, {"error": "r"}, {"error": "TypeError: x", "bug": True})
    assert msg is not None and "internal" in msg
