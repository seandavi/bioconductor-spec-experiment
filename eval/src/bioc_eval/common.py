"""Paths, JSONL io, case hashing, op vocabulary. Shared by every eval stage."""

from __future__ import annotations

import hashlib
import json
import tomllib
from pathlib import Path
from typing import Any, Iterable, Iterator

EVAL = Path(__file__).resolve().parents[2]
SPEC = EVAL / "spec"
CASES = EVAL / "cases"
BUILD = EVAL / "build"
RESULTS = EVAL / "results"
WORKTREES = Path(tomllib.loads((EVAL / "ports.toml").read_text())["WORKTREES"])

SOURCES = ("e1", "e1-redteam", "e2", "e4", "heldout")
POOL_FILES = ("e1", "e2", "e4", "heldout")


def canonical(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _ordered(node: Any) -> Any:
    """Hash view: sorted keys for unordered structural mappings, ordered pairs where order is semantic.

    Call `kwargs` (extra named arguments become mcols columns in call order) and `inputs` (evaluation order)
    are written as [name, value] pairs in their given order.
    """
    if isinstance(node, dict):
        out = {}
        for k, v in node.items():
            if k in ("kwargs", "inputs") and isinstance(v, dict):
                out[k] = [[kk, _ordered(vv)] for kk, vv in v.items()]
            else:
                out[k] = _ordered(v)
        return out
    if isinstance(node, list):
        return [_ordered(v) for v in node]
    return node


def case_hash(case: dict) -> str:
    """sha256 of canonical JSON of {inputs, call}; argument and input order are preserved (CASES.md, Hash)."""
    payload = _ordered({"inputs": case["inputs"], "call": case["call"]})
    return hashlib.sha256(canonical(payload).encode()).hexdigest()


def read_jsonl(path: Path) -> list[dict]:
    with path.open() as fh:
        return [json.loads(line) for line in fh if line.strip()]


def iter_jsonl(path: Path) -> Iterator[dict]:
    with path.open() as fh:
        for line in fh:
            if line.strip():
                yield json.loads(line)


def write_jsonl(path: Path, rows: Iterable[dict]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w") as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
            n += 1
    return n


def load_ops() -> dict[str, dict]:
    ops = json.loads((SPEC / "ops.json").read_text())["ops"]
    return {o["op"]: o for o in ops}


def load_families() -> dict[str, str]:
    return {name: o["family"] for name, o in load_ops().items()}


def ports() -> dict[str, dict]:
    return tomllib.loads((EVAL / "ports.toml").read_text())["ports"]


def load_pool() -> list[dict]:
    """All committed cases, sorted by id."""
    cases: list[dict] = []
    for stem in POOL_FILES:
        path = CASES / f"{stem}.jsonl"
        if path.exists():
            cases.extend(read_jsonl(path))
    cases.sort(key=lambda c: c["id"])
    return cases


def ops_in_call(call: dict) -> Iterator[str]:
    yield call["op"]
    for a in [*call["args"], *call["kwargs"].values()]:
        if isinstance(a, dict) and "$call" in a:
            yield from ops_in_call(a["$call"])
