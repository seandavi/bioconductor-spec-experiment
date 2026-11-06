"""Generate cases/grid_metadata.yaml: names / mcols / seqinfo propagation through every op.

    uv run python tests/differential/gen/grid_metadata.py

Each op runs on inputs carrying names, mcols (int, double, chr, lgl columns) and, for GRanges,
a seqinfo with lengths/circularity/genome. The canonical value includes all metadata, so any
op that drops or keeps it differently from R shows up.
"""

from __future__ import annotations

import random
from pathlib import Path

import yaml

rng = random.Random(7)


def mc(n: int) -> dict:
    return {"i": list(range(n)), "d": [round(rng.random(), 3) for _ in range(n)], "c": [rng.choice("xyz") for _ in range(n)],
            "l": [rng.random() < 0.5 for _ in range(n)]}


def ir(n: int, names: bool = True) -> dict:
    d = {"start": [rng.randint(1, 30) for _ in range(n)], "width": [rng.choice([0, 1, 3, 6]) for _ in range(n)], "mcols": mc(n)}
    if names:
        d["names"] = [f"r{i}" for i in range(n)]
    return {"IRanges": d}


SI = {"seqnames": ["chrA", "chrB", "chrC"], "seqlengths": [100, 200, None], "is_circular": [False, None, True], "genome": ["g1", "g1", None]}


def gr(n: int, names: bool = True, si: dict | None = SI) -> dict:
    d = {"seqnames": [rng.choice(["chrA", "chrB", "chrC"]) for _ in range(n)], "start": [rng.randint(1, 30) for _ in range(n)],
         "width": [rng.choice([0, 1, 3, 6]) for _ in range(n)], "strand": [rng.choice(["+", "-", "*"]) for _ in range(n)], "mcols": mc(n)}
    if si:
        d["seqinfo"] = si
    if names:
        d["names"] = [f"g{i}" for i in range(n)]
    return {"GRanges": d}


cases: list[dict] = []


def add(cid: str, inputs: dict, fn: str, args: list, kwargs: dict | None = None) -> None:
    call: dict = {"fn": fn, "args": args}
    if kwargs:
        call["kwargs"] = kwargs
    cases.append({"id": cid, "tags": ["metadata", fn], "inputs": inputs, "call": call})


OPS: list[tuple[str, list, dict]] = [
    ("shift", [3], {}), ("shift", [3], {"use.names": False}), ("narrow", [], {"start": 1, "end": 0}), ("narrow", [], {"start": 1, "width": 0}),
    ("resize", [2], {}), ("resize", [2], {"use.names": False}), ("flank", [2], {}), ("flank", [2], {"both": True, "use.names": False}),
    ("promoters", [1, 1], {}), ("sort", [], {}), ("sort", [], {"decreasing": True}), ("unique", [], {}), ("rev", [], {}),
    ("reduce", [], {}), ("reduce", [], {"with.revmap": True}), ("disjoin", [], {}), ("gaps", [], {"start": 1, "end": 40}), ("range", [], {}),
]
for kind, mk in (("ir", ir), ("gr", gr)):
    for k in range(2):
        x = mk(rng.choice([4, 7]))
        for fn, args, kw in OPS:
            if fn == "rev":
                continue
            suffix = "-" + "-".join(f"{a}{b}".replace(".", "") for a, b in kw.items()) if kw else ""
            add(f"meta-{kind}{k}-{fn}{suffix}", {"x": x}, fn, ["$x", *args], kw)
        n = len(x["IRanges" if kind == "ir" else "GRanges"]["start"])
        add(f"meta-{kind}{k}-subset-names", {"x": x}, "[", ["$x", [f"{'r' if kind == 'ir' else 'g'}{n - 1}", f"{'r' if kind == 'ir' else 'g'}0"]])
        add(f"meta-{kind}{k}-subset-logical", {"x": x}, "[", ["$x", [rng.random() < 0.5 for _ in range(n)]])
        add(f"meta-{kind}{k}-subset-negative", {"x": x}, "[", ["$x", [-1, -2]])
        add(f"meta-{kind}{k}-subset-repeat", {"x": x}, "[", ["$x", [1, 1, 2]])
        add(f"meta-{kind}{k}-names-set-null", {"x": x}, "names<-", ["$x"], {"value": {"NULL": []}})
        y = mk(rng.choice([3, 5]))
        add(f"meta-{kind}{k}-concat", {"x": x, "y": y}, "c", ["$x", "$y"])
        add(f"meta-{kind}{k}-concat-unnamed-y", {"x": x, "y": mk(3, names=False)}, "c", ["$x", "$y"])
        add(f"meta-{kind}{k}-concat-unnamed-x", {"x": mk(3, names=False), "y": y}, "c", ["$x", "$y"])
        add(f"meta-{kind}{k}-concat-self-twice", {"x": x}, "c", ["$x", "$x", "$x"])
        if kind == "gr":
            for fn in ("union", "intersect", "setdiff"):
                add(f"meta-gr{k}-{fn}", {"x": x, "y": y}, fn, ["$x", "$y"])
            add(f"meta-gr{k}-trim", {"x": x}, "trim", ["$x"])
            add(f"meta-gr{k}-subset-by-overlaps", {"x": x, "y": y}, "subsetByOverlaps", ["$x", "$y"])
            add(f"meta-gr{k}-subset-by-overlaps-invert", {"x": x, "y": y}, "subsetByOverlaps", ["$x", "$y"], {"invert": True})
            add(f"meta-gr{k}-ranges", {"x": x}, "ranges", ["$x"])
            add(f"meta-gr{k}-mcols", {"x": x}, "mcols", ["$x"])
            add(f"meta-gr{k}-seqinfo", {"x": x}, "seqinfo", ["$x"])
            add(f"meta-gr{k}-coverage", {"x": x}, "coverage", ["$x"])
        else:
            for fn in ("union", "intersect", "setdiff"):
                add(f"meta-ir{k}-{fn}", {"x": x, "y": y}, fn, ["$x", "$y"])
            add(f"meta-ir{k}-subset-by-overlaps", {"x": x, "y": y}, "subsetByOverlaps", ["$x", "$y"])
            add(f"meta-ir{k}-subset-by-overlaps-invert", {"x": x, "y": y}, "subsetByOverlaps", ["$x", "$y"], {"invert": True})
            add(f"meta-ir{k}-count", {"x": x, "y": y}, "countOverlaps", ["$x", "$y"])

out = Path(__file__).parent.parent / "cases" / "grid_metadata.yaml"
out.write_text("# GENERATED by gen/grid_metadata.py - edit the generator, not this file.\n" + yaml.safe_dump({"cases": cases}, default_flow_style=None, width=200))
print(f"wrote {out} ({len(cases)} cases)")
