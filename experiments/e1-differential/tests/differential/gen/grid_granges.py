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

for tag in ("r1", "r2", "r3"):
    x = gr(rng.choice([10, 20, 30]))
    y = gr(rng.choice([8, 20, 35]))
    xy = {"x": x, "y": y}
    for ign in (False, True):
        sfx = "-ign" if ign else ""
        kw0 = {"ignore.strand": True} if ign else {}
        cases.append({"id": f"grid-gr-{tag}-range{sfx}", "tags": ["grid", "range"], "inputs": {"x": x}, "call": {"fn": "range", "args": ["$x"], "kwargs": kw0}})
        for mg in (0, 1, 3):
            for de in (False, True):
                cases.append({"id": f"grid-gr-{tag}-reduce-g{mg}-{'drop' if de else 'keep'}{sfx}", "tags": ["grid", "reduce"], "inputs": {"x": x},
                              "call": {"fn": "reduce", "args": ["$x"], "kwargs": {"min.gapwidth": mg, "drop.empty.ranges": de, "with.revmap": True, **kw0}}})
        cases.append({"id": f"grid-gr-{tag}-disjoin{sfx}", "tags": ["grid", "disjoin"], "inputs": {"x": x}, "call": {"fn": "disjoin", "args": ["$x"], "kwargs": {"with.revmap": True, **kw0}}})
        cases.append({"id": f"grid-gr-{tag}-gaps{sfx}", "tags": ["grid", "gaps"], "inputs": {"x": x}, "call": {"fn": "gaps", "args": ["$x"], "kwargs": kw0}})
        cases.append({"id": f"grid-gr-{tag}-gaps-window{sfx}", "tags": ["grid", "gaps"], "inputs": {"x": x}, "call": {"fn": "gaps", "args": ["$x"], "kwargs": {"start": 3, "end": 30, **kw0}}})
        for fn in ("union", "intersect", "setdiff"):
            cases.append({"id": f"grid-gr-{tag}-{fn}{sfx}", "tags": ["grid", fn], "inputs": xy, "call": {"fn": fn, "args": ["$x", "$y"], "kwargs": kw0}})
        for fix in ("start", "end", "center"):
            cases.append({"id": f"grid-gr-{tag}-resize-{fix}{sfx}", "tags": ["grid", "resize"], "inputs": {"x": x}, "call": {"fn": "resize", "args": ["$x", rng.randint(0, 9)], "kwargs": {"fix": fix, **kw0}}})
        cases.append({"id": f"grid-gr-{tag}-flank{sfx}", "tags": ["grid", "flank"], "inputs": {"x": x},
                      "call": {"fn": "flank", "args": ["$x", rng.randint(-5, 6)], "kwargs": {"start": rng.random() < 0.5, "both": rng.random() < 0.3, **kw0}}})
    cases.append({"id": f"grid-gr-{tag}-promoters", "tags": ["grid", "promoters"], "inputs": {"x": x}, "call": {"fn": "promoters", "args": ["$x", rng.randint(0, 6), rng.randint(0, 6)]}})
    cases.append({"id": f"grid-gr-{tag}-trim-promoters", "tags": ["grid", "trim"], "inputs": {"x": x}, "call": {"fn": "trim", "args": [{"call": {"fn": "promoters", "args": ["$x", 40, 5]}}]}})
    cases.append({"id": f"grid-gr-{tag}-shift", "tags": ["grid", "shift"], "inputs": {"x": x}, "call": {"fn": "shift", "args": ["$x", rng.randint(-10, 10)]}})

out = Path(__file__).parent.parent / "cases" / "grid_granges.yaml"
out.write_text("# GENERATED by gen/grid_granges.py - edit the generator, not this file.\n" + yaml.safe_dump({"cases": cases}, default_flow_style=None, width=200))
print(f"wrote {out} ({len(cases)} cases)")
