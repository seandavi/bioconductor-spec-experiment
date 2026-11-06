"""R-result cache lifecycle (finding 1): one fingerprint is written, compared, and tracks script changes."""

import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "differential"))
import harness  # noqa: E402

CASE = harness.Case("c1", "unit_cache", {"id": "c1", "inputs": {}, "call": {"fn": "identity", "args": [1]}})
RV = "R 9.9.9; IRanges 1; GenomicRanges 2; S4Vectors 3"


@pytest.fixture
def cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    state: dict[str, Any] = {"r_calls": 0, "r": True}
    monkeypatch.setattr(harness, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(harness, "r_available", lambda: state["r"])
    monkeypatch.setattr(harness, "r_versions", lambda: RV)
    monkeypatch.setattr(harness, "_fresh", {})

    def fake_run_r(cases: list[harness.Case]) -> tuple[dict[str, Any], str]:
        state["r_calls"] += 1
        return {c.id: {"value": 1, "warnings": []} for c in cases}, RV

    monkeypatch.setattr(harness, "run_r", fake_run_r)
    return state


def reload() -> None:
    harness._fresh.clear()  # noqa: SLF001 - simulate a new pytest session


def test_refreshed_cache_is_a_hit_in_the_next_session(cache: dict[str, Any]) -> None:
    assert harness.r_result(CASE, [CASE]) is not None and cache["r_calls"] == 1
    reload()
    assert harness.r_result(CASE, [CASE]) is not None and cache["r_calls"] == 1  # no rerun


def test_persisted_fingerprint_is_the_compared_one(cache: dict[str, Any]) -> None:
    harness.r_result(CASE, [CASE])
    assert harness._load_cache("unit_cache")["versions"] == harness.fingerprint(RV)  # noqa: SLF001


def test_script_change_invalidates_online_and_offline(cache: dict[str, Any], monkeypatch: pytest.MonkeyPatch) -> None:
    harness.r_result(CASE, [CASE])
    monkeypatch.setattr(harness, "script_digest", lambda: "changed")
    reload()
    assert harness.r_result(CASE, [CASE]) is not None and cache["r_calls"] == 2  # online: rerun
    monkeypatch.setattr(harness, "script_digest", lambda: "changed-again")
    cache["r"] = False
    reload()
    assert harness.r_result(CASE, [CASE]) is None  # offline: stale, skipped rather than trusted


def test_version_change_invalidates(cache: dict[str, Any], monkeypatch: pytest.MonkeyPatch) -> None:
    harness.r_result(CASE, [CASE])
    monkeypatch.setattr(harness, "r_versions", lambda: RV.replace("9.9.9", "10.0.0"))
    reload()
    harness.r_result(CASE, [CASE])
    assert cache["r_calls"] == 2


def test_offline_trusts_cache_with_matching_scripts(cache: dict[str, Any]) -> None:
    harness.r_result(CASE, [CASE])
    cache["r"] = False
    reload()
    assert harness.r_result(CASE, [CASE]) is not None and cache["r_calls"] == 1
