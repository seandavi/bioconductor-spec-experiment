"""Generate cases/grid_granges.yaml: random small GRanges pairs through overlap/nearest ops.

    uv run python tests/differential/gen/grid_granges.py

Deterministic (fixed seed). Coordinates are tiny so ties, adjacency, empty ranges and
strand combinations are frequent. Sizes cover nq < ns, nq == ns and nq > ns.
"""

from __future__ import annotations

import random
from pathlib import Path

import yaml

rng = random.Random(20261105)
SEQINFO = {"seqnames": ["chr1", "chr2", "chr3"], "seqlengths": [1000, 1000, 1000]}


def gr(n: int, seqs: tuple[str, ...] = ("chr1", "chr2"), span: int = 25) -> dict:
    return {"GRanges": {
        "seqnames": [rng.choice(seqs) for _ in range(n)],
        "start": [rng.randint(1, span) for _ in range(n)],
        "width": [rng.choice([0, 1, 1, 2, 3, 5, 8]) for _ in range(n)],
        "strand": [rng.choice(["+", "-", "*"]) for _ in range(n)],
        "seqinfo": SEQINFO,
    }}


cases: list[dict] = []
SIZES = {"qshort": (6, 30), "equal": (20, 20), "qlong": (35, 12)}
for tag, (nq, ns) in SIZES.items():
    q, s = gr(nq), gr(ns)
    for typ in ("any", "start", "end", "within", "equal"):
        for ign in (False, True):
            kw: dict = {"type": typ}
            if ign:
                kw["ignore.strand"] = True
            cases.append({"id": f"grid-gr-fo-{tag}-{typ}{'-ignstrand' if ign else ''}", "tags": ["grid", "findOverlaps", typ],
                          "inputs": {"q": q, "s": s}, "call": {"fn": "findOverlaps", "args": ["$q", "$s"], "kwargs": kw}})
    for sel in ("first", "last", "arbitrary"):
        cases.append({"id": f"grid-gr-fo-{tag}-select-{sel}", "tags": ["grid", "findOverlaps", "select"],
                      "inputs": {"q": q, "s": s}, "call": {"fn": "findOverlaps", "args": ["$q", "$s"], "kwargs": {"select": sel}}})
    cases.append({"id": f"grid-gr-fo-{tag}-maxgap", "tags": ["grid", "findOverlaps"],
                  "inputs": {"q": q, "s": s}, "call": {"fn": "findOverlaps", "args": ["$q", "$s"], "kwargs": {"maxgap": 2, "type": "within"}}})
    cases.append({"id": f"grid-gr-count-{tag}", "tags": ["grid", "countOverlaps"],
                  "inputs": {"q": q, "s": s}, "call": {"fn": "countOverlaps", "args": ["$q", "$s"], "kwargs": {"minoverlap": 2}}})
    for fn, sels in (("nearest", ("arbitrary", "all")), ("precede", ("first", "all")), ("follow", ("last", "all"))):
        for sel in sels:
            for ign in (False, True):
                kw = {"select": sel}
                if ign:
                    kw["ignore.strand"] = True
                cases.append({"id": f"grid-gr-{fn}-{tag}-{sel}{'-ignstrand' if ign else ''}", "tags": ["grid", fn],
                              "inputs": {"q": q, "s": s}, "call": {"fn": fn, "args": ["$q", "$s"], "kwargs": kw}})
    d1, d2 = gr(nq), gr(nq)
    cases.append({"id": f"grid-gr-distance-{tag}", "tags": ["grid", "distance"], "inputs": {"x": d1, "y": d2},
                  "call": {"fn": "distance", "args": ["$x", "$y"]}})
    cases.append({"id": f"grid-gr-distance-ignstrand-{tag}", "tags": ["grid", "distance"], "inputs": {"x": d1, "y": d2},
                  "call": {"fn": "distance", "args": ["$x", "$y"], "kwargs": {"ignore.strand": True}}})
for tag, n in (("a", 12), ("b", 25)):
    x = gr(n)
    for fn in ("nearest", "precede", "follow"):
        for sel in ({"nearest": ("arbitrary", "all")}.get(fn) or ({"precede": ("first", "all"), "follow": ("last", "all")}[fn])):
            cases.append({"id": f"grid-gr-{fn}-self-{tag}-{sel}", "tags": ["grid", fn, "self"], "inputs": {"x": x},
                          "call": {"fn": fn, "args": ["$x"], "kwargs": {"select": sel}}})
    for kw_name, kw in (("plain", {}), ("dropself", {"drop.self": True}), ("dropboth", {"drop.self": True, "drop.redundant": True})):
        cases.append({"id": f"grid-gr-fo-self-{tag}-{kw_name}", "tags": ["grid", "findOverlaps", "self"], "inputs": {"x": x},
                      "call": {"fn": "findOverlaps", "args": ["$x"], "kwargs": kw}})
    cases.append({"id": f"grid-gr-order-{tag}", "tags": ["grid", "order"], "inputs": {"x": x}, "call": {"fn": "order", "args": ["$x"]}})
    cases.append({"id": f"grid-gr-order-decreasing-{tag}", "tags": ["grid", "order"], "inputs": {"x": x}, "call": {"fn": "order", "args": ["$x"], "kwargs": {"decreasing": True}}})
    for tm in ("min", "max", "first", "last", "average"):
        cases.append({"id": f"grid-gr-rank-{tm}-{tag}", "tags": ["grid", "rank"], "inputs": {"x": x}, "call": {"fn": "rank", "args": ["$x"], "kwargs": {"ties.method": tm}}})
    cases.append({"id": f"grid-gr-duplicated-{tag}", "tags": ["grid", "duplicated"], "inputs": {"x": x}, "call": {"fn": "duplicated", "args": ["$x"]}})
    cases.append({"id": f"grid-gr-duplicated-fromlast-{tag}", "tags": ["grid", "duplicated"], "inputs": {"x": x}, "call": {"fn": "duplicated", "args": ["$x"], "kwargs": {"fromLast": True}}})
    cases.append({"id": f"grid-gr-unique-{tag}", "tags": ["grid", "unique"], "inputs": {"x": x}, "call": {"fn": "unique", "args": ["$x"]}})
    cases.append({"id": f"grid-gr-sort-{tag}", "tags": ["grid", "sort"], "inputs": {"x": x}, "call": {"fn": "sort", "args": ["$x"]}})
    cases.append({"id": f"grid-gr-sort-ignstrand-{tag}", "tags": ["grid", "sort"], "inputs": {"x": x}, "call": {"fn": "sort", "args": ["$x"], "kwargs": {"ignore.strand": True}}})
    cases.append({"id": f"grid-gr-coverage-{tag}", "tags": ["grid", "coverage"], "inputs": {"x": x}, "call": {"fn": "coverage", "args": ["$x"]}})

out = Path(__file__).parent.parent / "cases" / "grid_granges.yaml"
out.write_text("# GENERATED by gen/grid_granges.py - edit the generator, not this file.\n" + yaml.safe_dump({"cases": cases}, default_flow_style=None, width=200))
print(f"wrote {out} ({len(cases)} cases)")
