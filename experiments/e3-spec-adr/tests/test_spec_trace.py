"""Spec/test traceability: every live clause is cited by a test, every cited id exists.

Clauses are ``**ID.n**`` at the start of a bullet-free line in ``spec/*.md``; a
retired clause carries ``Retired:`` on its line. Tests cite ids with
``@pytest.mark.spec("ID.n", ...)``.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLAUSE = re.compile(r"^\*\*([A-Z]+\.\d+)\*\*(.*)$", re.M)

# Areas whose implementation has not landed yet (shrinks each stage; empty at close-out).
PENDING_AREAS: set[str] = {"CONV", "INTRA", "INTER", "SET", "COV", "OVL", "NEAR", "ORD", "SEQI", "GR"}


def spec_clauses() -> dict[str, bool]:
    """Clause id -> is_live."""
    out: dict[str, bool] = {}
    for path in sorted((ROOT / "spec").glob("*.md")):
        for m in CLAUSE.finditer(path.read_text()):
            assert m.group(1) not in out, f"duplicate clause id {m.group(1)}"
            out[m.group(1)] = "Retired:" not in m.group(2)
    return out


def cited_ids() -> set[str]:
    ids: set[str] = set()
    for path in sorted((ROOT / "tests").glob("test_*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call) and _is_spec_mark(node.func):
                ids.update(a.value for a in node.args if isinstance(a, ast.Constant))
    return ids


def _is_spec_mark(func: ast.expr) -> bool:
    return isinstance(func, ast.Attribute) and func.attr == "spec" and (
        isinstance(func.value, ast.Attribute) and func.value.attr == "mark"
    )


def test_every_live_clause_is_cited() -> None:
    clauses = spec_clauses()
    cited = cited_ids()
    missing = sorted(
        c for c, live in clauses.items()
        if live and c not in cited and c.split(".")[0] not in PENDING_AREAS
    )
    assert not missing, f"clauses without a test: {missing}"


def test_every_cited_id_exists() -> None:
    unknown = sorted(cited_ids() - set(spec_clauses()))
    assert not unknown, f"tests cite unknown clauses: {unknown}"
