"""Join resolved cases, R expected, raw port results and scored outcomes (stdlib only)."""
from __future__ import annotations

import json
from pathlib import Path

EVAL = Path(__file__).resolve().parent.parent
PORTS = ("e0", "e0b", "e1", "e2", "e3", "e4", "biocpy")
ADJUDICATED = ("mismatch", "port_error_unexpected")


def _jsonl(path: Path) -> dict[str, dict]:
    with path.open() as fh:
        return {r["id"]: r for r in map(json.loads, fh)}


def load() -> tuple[dict[str, dict], dict[str, dict], dict[str, dict[str, dict]], dict[str, dict[str, dict]]]:
    resolved = _jsonl(EVAL / "build" / "resolved.jsonl")
    expected = _jsonl(EVAL / "build" / "expected.jsonl")
    raw = {p: _jsonl(EVAL / "build" / "results" / f"{p}.jsonl") for p in PORTS}
    scored = {p: _jsonl(EVAL / "results" / "scored" / f"{p}.jsonl") for p in PORTS}
    return resolved, expected, raw, scored
