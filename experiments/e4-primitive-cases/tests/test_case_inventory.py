"""Consistency of clauses <-> cases <-> goldens <-> ADRs (PLAN 4.5)."""

from __future__ import annotations

import re

import pytest

from conformance.loader import ADR_DIR, ROOT, golden_ids, load_cases

CASES = load_cases()

# TARGET.md operations, each mapped to a clause-id substring that must exist.
TARGET_OPS = {
    "shift": "-SHIFT-", "narrow": "-NARROW-", "resize": "-RESIZE-", "flank": "-FLANK-",
    "promoters": "-PROMOTERS-", "range": "-RANGE-", "reduce": "-REDUCE-", "gaps": "-GAPS-",
    "disjoin": "-DISJOIN-", "union": "-UNION-", "intersect": "-INTERSECT-",
    "setdiff": "-SETDIFF-", "findOverlaps": "-FO-", "countOverlaps": "-COUNT-",
    "overlapsAny": "-ANY-", "subsetByOverlaps": "-SUBSET-", "nearest": "-NEAREST-",
    "precede": "-PRECEDE-", "follow": "-FOLLOW-", "distance": "-DISTANCE-",
    "coverage": "-COV-", "order": "-ORDER-", "sort": "-SORT-", "rank": "-RANK-",
    "duplicated": "-DUP-", "unique": "-UNIQUE-", "trim": "-TRIM-",
    "seqinfo merge": "-MERGE-", "ignore.strand": "-IGNSTRAND-",
}


def clause_ids() -> set[str]:
    text = (ROOT / "clauses.md").read_text()
    return set(re.findall(r"^\| ((?:IR|GR|SI)-[A-Z0-9-]+) \|", text, flags=re.M))


def adr_ids() -> set[str]:
    return {f"ADR-{p.name[:4]}" for p in ADR_DIR.glob("[0-9]*.md")}


def test_clause_ids_unique() -> None:
    ids = re.findall(r"^\| ((?:IR|GR|SI)-[A-Z0-9-]+) \|", (ROOT / "clauses.md").read_text(), flags=re.M)
    assert len(ids) == len(set(ids))


def test_case_ids_unique_and_prefixed() -> None:
    ids = [c.id for c in CASES]
    assert len(ids) == len(set(ids))
    for c in CASES:
        assert c.id.startswith(("ir-", "gr-")), c.id
        assert c.id.startswith(c.source.parent.name[:2] + "-"), c.id


def test_every_case_names_an_existing_clause() -> None:
    known = clause_ids()
    assert [c.id for c in CASES if c.spec["clause"] not in known] == []


def test_every_clause_has_a_case() -> None:
    used = {c.spec["clause"] for c in CASES}
    assert sorted(clause_ids() - used) == []


def test_goldens_match_cases() -> None:
    assert {c.id for c in CASES} == golden_ids()


def test_goldens_have_no_expected_values_in_cases() -> None:
    forbidden = {"expect", "expected", "golden", "result"}
    for c in CASES:
        assert not forbidden & set(c.spec), c.id


@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_non_exact_modes_and_overrides_cite_an_adr(case) -> None:  # type: ignore[no-untyped-def]
    cmp = case.spec.get("compare") or {}
    non_exact = bool(cmp.get("rows")) or bool(cmp.get("fields"))
    if non_exact or case.spec.get("python"):
        adr = case.spec.get("adr")
        assert adr, f"{case.id}: non-exact compare or python override needs `adr:`"
        assert adr in adr_ids(), f"{case.id}: {adr} has no file in adr/"


def test_target_operations_covered_by_clauses() -> None:
    ids = " ".join(clause_ids())
    missing = [op for op, frag in TARGET_OPS.items() if frag not in ids]
    assert missing == [] or not (ROOT / "STAGE_DONE").exists(), missing
