"""Consistency of clauses <-> cases <-> goldens <-> ADRs (PLAN 4.5)."""

from __future__ import annotations

import re

import pytest

from conformance.loader import ADR_DIR, ROOT, golden_ids, load_cases

CASES = load_cases()

# TARGET.md operations: each must be named in a clause id or statement.
TARGET_OPS = {
    "shift": "IR-SHIFT-", "narrow": "IR-NARROW-", "resize": "IR-RESIZE-", "flank": "IR-FLANK-",
    "promoters": "IR-PROMOTERS-", "range": "IR-RANGE-", "reduce": "IR-REDUCE-", "gaps": "IR-GAPS-",
    "disjoin": "IR-DISJOIN-", "union": "IR-UNION-", "intersect": "IR-INTERSECT-",
    "setdiff": "IR-SETDIFF-", "punion": "IR-PUNION-", "pintersect": "IR-PINTERSECT-",
    "psetdiff": "IR-PSETDIFF-", "findOverlaps": "IR-FO-", "countOverlaps": "IR-COUNT-",
    "overlapsAny": "IR-OVANY-", "subsetByOverlaps": "IR-SUBSET-", "nearest": "IR-NEAREST-",
    "precede": "IR-PRECEDE-", "follow": "IR-FOLLOW-", "distance": "IR-DISTANCE-",
    "coverage": "IR-COV-", "order": "IR-ORDER-", "sort": "IR-SORT-", "rank": "IR-RANK-",
    "duplicated": "IR-DUP-", "unique": "IR-UNIQUE-", "subsetting": "IR-SUB-",
    "concatenation": "IR-CAT-", "accessors": "IR-ACC-", "construction": "IR-CON-",
    "mcols": "IR-CON-MCOLS-", "GRanges construction": "GR-CON-", "GRanges accessors": "GR-ACC-",
    "Seqinfo": "SI-CON-", "seqinfo merge": "GR-CAT-1", "trim": "GR-TRIM-",
    "GRanges coverage": "GR-COV-", "GRanges ordering": "GR-ORDER-", "ignore.strand": "ignore.strand",
}

# Operations that must be exercised by at least one case per class (call op names).
CASE_OPS_BOTH = {
    "shift", "narrow", "resize", "flank", "promoters", "range", "reduce", "gaps", "disjoin",
    "union", "intersect", "setdiff", "punion", "pintersect", "psetdiff", "findOverlaps",
    "countOverlaps", "overlapsAny", "subsetByOverlaps", "nearest", "precede", "follow",
    "distance", "coverage", "order", "sort", "rank", "duplicated", "unique", "[", "c",
    "start", "end", "width", "names", "length",
}
CASE_OPS_GRANGES_ONLY = {"seqnames", "strand", "ranges", "seqlevels", "seqlengths", "isCircular", "genome", "seqinfo", "trim"}


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
    text = (ROOT / "clauses.md").read_text()
    missing = [op for op, frag in TARGET_OPS.items() if frag not in text]
    assert missing == []


@pytest.mark.parametrize("prefix,ops", [("ir-", CASE_OPS_BOTH), ("gr-", CASE_OPS_BOTH | CASE_OPS_GRANGES_ONLY)])
def test_every_target_operation_has_cases(prefix: str, ops: set[str]) -> None:
    used = {c.spec["call"]["op"] for c in CASES if c.id.startswith(prefix)}
    assert sorted(ops - used) == []
