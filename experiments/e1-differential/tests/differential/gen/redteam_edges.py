"""Build empty/zero-width, negative-coordinate, and ordering adversaries."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

OUT = Path(__file__).parents[1] / "cases" / "redteam_edges.yaml"
cases: list[dict[str, Any]] = []


def add(cid: str, inputs: dict[str, Any], fn: str, args: list[Any] | None = None,
        **kwargs: Any) -> None:
    cases.append({"id": f"rt-edge-{cid}", "tags": ["redteam", "edges", fn], "inputs": inputs,
                  "call": {"fn": fn, "args": args if args is not None else ["$x"], "kwargs": kwargs}})


def raw(cid: str, r: str, py: str) -> None:
    cases.append({"id": f"rt-edge-{cid}", "tags": ["redteam", "edges", "construct"], "raw": {"r": r, "py": py}})


raw("ir-default-empty", "IRanges()", "IRanges()")
raw("gr-default-empty", "GRanges()", "GRanges()")
empty_ir: dict[str, Any] = {"IRanges": {"start": [], "width": [], "names": {"chr": []}, "mcols": {"score": {"int": []}}}}
empty_gr: dict[str, Any] = {"GRanges": {"seqnames": {"chr": []}, "start": [], "width": [], "strand": {"chr": []},
    "names": {"chr": []}, "mcols": {"score": {"int": []}},
    "seqinfo": {"seqnames": ["chrZ", "chrA"], "seqlengths": [0, 4], "is_circular": [False, False]}}}
for tag, x in (("ir", empty_ir), ("gr", empty_gr)):
    for fn, kwargs in (("reduce", {"with.revmap": True}), ("disjoin", {"with.revmap": True}),
                       ("range", {}), ("coverage", {}), ("sort", {}), ("unique", {}), ("nearest", {"select": "all"})):
        add(f"{tag}-empty-{fn}", {"x": x}, fn, **kwargs)
    add(f"{tag}-empty-union", {"x": x, "y": x}, "union", ["$x", "$y"])

x = {"GRanges": {"seqnames": ["chrZ", "chrZ", "chrZ", "chrA", "chrA"],
                 "start": [-2, 1, 3, -2, 4], "width": [2, 0, 0, 3, 0],
                 "strand": ["+", "-", "*", "+", "*"], "names": ["a", "b", "c", "d", "e"],
                 "mcols": {"score": [1, 2, 3, 4, 5]},
                 "seqinfo": {"seqnames": ["chrZ", "chrA", "unused"], "seqlengths": [0, 3, 2],
                             "is_circular": [False, False, False]}}}
add("zero-seqlength-trim", {"x": x}, "trim")
add("zero-seqlength-gaps", {"x": x}, "gaps")
add("zero-seqlength-coverage", {"x": x}, "coverage")

# Ties differ only in metadata and names. Seqlevels deliberately oppose lexical
# order; '*' is an ordering level, not a wildcard for duplication.
x = {"GRanges": {"seqnames": ["chrA", "chrZ", "chrZ", "chrZ", "chrZ", "chrA"],
                 "start": [5, 5, 5, 5, 5, 5], "width": [0, 0, 0, 0, 0, 0],
                 "strand": ["*", "*", "-", "+", "+", "*"], "names": ["a", "b", "c", "d", "e", "f"],
                 "mcols": {"score": [6, 5, 4, 3, 2, 1]}, "seqinfo": {"seqnames": ["chrZ", "chrA", "unused"]}}}
for fn, suffix, kwargs in (
    ("order", "ascending", {}), ("order", "descending", {"decreasing": True}),
    ("sort", "descending", {"decreasing": True}), ("sort", "ignore-strand", {"ignore.strand": True}),
    ("rank", "last", {"ties.method": "last"}),
    ("duplicated", "from-last", {"fromLast": True}), ("unique", "from-last", {"fromLast": True}),
):
    add(f"ordering-{fn}-{suffix}", {"x": x}, fn, **kwargs)
add("ordering-subset-empty-name", {"x": {"IRanges": {"start": [3, 1, 2], "width": [0, 1, 2],
    "names": ["", "dup", "dup"], "mcols": {"score": [1, 2, 3]}}}}, "[", ["$x", ["", "dup", "dup"]])

for tag, kind in (("ir", "IRanges"), ("gr", "GRanges")):
    q: dict[str, Any] = {"start": [-5, -1, 1, 3], "width": [0, 2, 0, 1]}
    s: dict[str, Any] = {"start": [-7, -3, 0, 5], "width": [1, 1, 1, 0]}
    if tag == "gr":
        q.update(seqnames=["chr2", "chr1", "chr2", "chr1"], strand=["+", "*", "-", "+"],
                 seqinfo={"seqnames": ["chr1", "chr2"]})
        s.update(seqnames=["chr1", "chr2", "chr1", "chr2"], strand=["+", "-", "*", "+"],
                 seqinfo={"seqnames": ["chr2", "chr1"]})
    for fn in ("nearest", "precede", "follow"):
        add(f"{tag}-negative-{fn}", {"x": {kind: q}, "y": {kind: s}}, fn, ["$x", "$y"], select="all")
    weight_x: dict[str, Any] = {"start": [1, 2, 5], "width": [3, 0, 2]}
    if tag == "gr":
        weight_x.update(seqnames=["chr1"])
    add(f"{tag}-coverage-typed-na", {"x": {kind: weight_x}}, "coverage", weight={"int": [1, None, 3]})
    add(f"{tag}-coverage-overlap-na", {"x": {kind: weight_x}}, "coverage", weight={"num": [None, 2, 3]})

# Compare unknown circularity with explicit linear/circular sequences of known length.
for circular, suffix in ((None, "unknown"), (False, "linear"), (True, "circular")):
    add(f"trim-circularity-{suffix}", {"x": {"GRanges": {"seqnames": ["chr1"], "start": [-2], "width": [8],
        "seqinfo": {"seqnames": ["chr1"], "seqlengths": [3], "is_circular": [circular]}}}}, "trim")

# Empty ranges alone still define extent for some R operations.
x = {"IRanges": {"start": [8, -2, 8, 4], "width": [0, 0, 0, 0], "names": ["a", "b", "c", "d"],
                 "mcols": {"score": [1, 2, 3, 4]}}}
for fn, kwargs in (("range", {}), ("gaps", {"start": -5, "end": 10}), ("reduce", {"with.revmap": True}),
                   ("disjoin", {"with.revmap": True})):
    add(f"all-zero-{fn}", {"x": x}, fn, **kwargs)

old = yaml.safe_load(OUT.read_text()) if OUT.exists() else {}
reasons = {c["id"]: c["xfail"] for c in old.get("cases", []) if "xfail" in c}
for case in cases:
    case["check_warnings"] = True
    if case["id"] in reasons:
        case["xfail"] = reasons[case["id"]]
OUT.write_text("# Generated by gen/redteam_edges.py; reviewed xfail reasons are retained.\n" +
               yaml.safe_dump({"cases": cases}, sort_keys=False, width=110))
print(f"wrote {len(cases)} cases")
