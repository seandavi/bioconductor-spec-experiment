"""Generate cases/grid_bad_args.yaml: invalid / odd arguments to every IRanges and GRanges op.

    uv run python tests/differential/gen/grid_bad_args.py

Each case passes R and Python the same malformed argument; both must agree on error vs value
(and on the value when there is one).
"""

from __future__ import annotations

from pathlib import Path

import yaml

X = {"IRanges": {"start": [1, 5, 9], "width": [3, 0, 2], "names": ["a", "b", "c"]}}
Y = {"IRanges": {"start": [2, 6], "width": [4, 1]}}
G = {"GRanges": {"seqnames": ["chr1", "chr1", "chr2"], "start": [1, 5, 9], "width": [3, 0, 2], "strand": ["+", "-", "*"],
                 "seqinfo": {"seqnames": ["chr1", "chr2"], "seqlengths": [100, 100]}}}
H = {"GRanges": {"seqnames": ["chr1", "chr2"], "start": [2, 6], "width": [4, 1], "strand": ["+", "*"],
                 "seqinfo": {"seqnames": ["chr1", "chr2"], "seqlengths": [100, 100]}}}
NA, STR, NUM2 = [None], ["a"], [1, 2]
VALS: dict[str, object] = {"NA": NA, "str": STR, "len2": NUM2, "neg": [-1], "dbl": {"num": [1.5]}, "lgl": [True], "empty": {"int": []}}
# 2147483647 is left out on purpose: R's integer arithmetic then wraps silently (invalid objects,
# missing hits); Python raises or answers exactly (ADR-0005).

cases: list[dict] = []


def add(cid: str, inputs: dict, fn: str, args: list, kwargs: dict | None = None) -> None:
    call: dict = {"fn": fn, "args": args}
    if kwargs:
        call["kwargs"] = kwargs
    cases.append({"id": f"bad-{cid}", "tags": ["bad-args", fn], "inputs": inputs, "call": call})


for tag, obj in (("ir", X), ("gr", G)):
    inp = {"x": obj, "y": Y if tag == "ir" else H}
    for vname, v in VALS.items():
        add(f"{tag}-shift-{vname}", inp, "shift", ["$x", v])
        add(f"{tag}-resize-width-{vname}", inp, "resize", ["$x", v])
        add(f"{tag}-flank-width-{vname}", inp, "flank", ["$x", v])
        add(f"{tag}-flank-start-{vname}", inp, "flank", ["$x", 2], {"start": v})
        add(f"{tag}-flank-both-{vname}", inp, "flank", ["$x", 2], {"both": v})
        add(f"{tag}-promoters-up-{vname}", inp, "promoters", ["$x", v, 1])
        add(f"{tag}-promoters-down-{vname}", inp, "promoters", ["$x", 1, v])
        add(f"{tag}-narrow-start-{vname}", inp, "narrow", ["$x"], {"start": v})
        add(f"{tag}-narrow-end-{vname}", inp, "narrow", ["$x"], {"end": v})
        add(f"{tag}-narrow-width-{vname}", inp, "narrow", ["$x"], {"start": 1, "width": v})
        add(f"{tag}-reduce-mingap-{vname}", inp, "reduce", ["$x"], {"min.gapwidth": v})
        add(f"{tag}-reduce-drop-{vname}", inp, "reduce", ["$x"], {"drop.empty.ranges": v})
        add(f"{tag}-reduce-revmap-{vname}", inp, "reduce", ["$x"], {"with.revmap": v})
        add(f"{tag}-gaps-start-{vname}", inp, "gaps", ["$x"], {"start": v})
        add(f"{tag}-gaps-end-{vname}", inp, "gaps", ["$x"], {"end": v})
        add(f"{tag}-disjoin-revmap-{vname}", inp, "disjoin", ["$x"], {"with.revmap": v})
        add(f"{tag}-fo-maxgap-{vname}", inp, "findOverlaps", ["$x", "$y"], {"maxgap": v})
        add(f"{tag}-fo-minoverlap-{vname}", inp, "findOverlaps", ["$x", "$y"], {"minoverlap": v, "type": "within"})
        add(f"{tag}-fo-type-{vname}", inp, "findOverlaps", ["$x", "$y"], {"type": v})
        add(f"{tag}-fo-select-{vname}", inp, "findOverlaps", ["$x", "$y"], {"select": v})
        add(f"{tag}-count-maxgap-{vname}", inp, "countOverlaps", ["$x", "$y"], {"maxgap": v})
        add(f"{tag}-nearest-select-{vname}", inp, "nearest", ["$x", "$y"], {"select": v})
        add(f"{tag}-precede-select-{vname}", inp, "precede", ["$x", "$y"], {"select": v})
        add(f"{tag}-follow-select-{vname}", inp, "follow", ["$x", "$y"], {"select": v})
        add(f"{tag}-coverage-shift-{vname}", inp, "coverage", ["$x"], {"shift": v})
        add(f"{tag}-coverage-width-{vname}", inp, "coverage", ["$x"], {"width": v})
        add(f"{tag}-coverage-weight-{vname}", inp, "coverage", ["$x"], {"weight": v})
        add(f"{tag}-order-decreasing-{vname}", inp, "order", ["$x"], {"decreasing": v})
        add(f"{tag}-sort-decreasing-{vname}", inp, "sort", ["$x"], {"decreasing": v})
        add(f"{tag}-rank-ties-{vname}", inp, "rank", ["$x"], {"ties.method": v})
        add(f"{tag}-duplicated-fromLast-{vname}", inp, "duplicated", ["$x"], {"fromLast": v})
        add(f"{tag}-unique-fromLast-{vname}", inp, "unique", ["$x"], {"fromLast": v})
        add(f"{tag}-subset-{vname}", inp, "[", ["$x", v])
        if tag == "ir":
            add(f"ir-punion-fillgap-{vname}", inp, "punion", ["$x", "$x"], {"fill.gap": v})
            add(f"ir-pintersect-resolve-{vname}", inp, "pintersect", ["$x", "$x"], {"resolve.empty": v})
        else:
            add(f"gr-resize-ignore-{vname}", inp, "resize", ["$x", 2], {"ignore.strand": v})
            add(f"gr-reduce-ignore-{vname}", inp, "reduce", ["$x"], {"ignore.strand": v})
            add(f"gr-fo-ignore-{vname}", inp, "findOverlaps", ["$x", "$y"], {"ignore.strand": v})
            add(f"gr-gaps-ignore-{vname}", inp, "gaps", ["$x"], {"ignore.strand": v})
            add(f"gr-union-ignore-{vname}", inp, "union", ["$x", "$y"], {"ignore.strand": v})
            add(f"gr-nearest-ignore-{vname}", inp, "nearest", ["$x", "$y"], {"ignore.strand": v})
            add(f"gr-distance-ignore-{vname}", inp, "distance", ["$x", "$x"], {"ignore.strand": v})
            add(f"gr-punion-fillgap-{vname}", inp, "punion", ["$x", "$x"], {"fill.gap": v})
            add(f"gr-pintersect-drop-{vname}", inp, "pintersect", ["$x", "$x"], {"drop.nohit.ranges": v})

out = Path(__file__).parent.parent / "cases" / "grid_bad_args.yaml"
out.write_text("# GENERATED by gen/grid_bad_args.py - edit the generator, not this file.\n" + yaml.safe_dump({"cases": cases}, default_flow_style=None, width=200))
print(f"wrote {out} ({len(cases)} cases)")
