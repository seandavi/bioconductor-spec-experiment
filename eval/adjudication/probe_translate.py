"""R probe: are nearest/precede/follow/distance and overlap results invariant under translation?

Those operations depend only on relative positions within a sequence, so adding K to every start
and end of every range input must not change the result (for index- and count-valued results).
Every such case without a circular sequence of known length is re-run in R with K = 1000 added; a different result means
R's expected value is an artefact of where the ranges sit (e.g. the sentinel number line GRanges
precede/follow/nearest build, which coordinates <= 0 cross).
Writes probe_translate.jsonl: {id, translation_dependent (None: shifted input overflows), r_shifted}.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

from load import EVAL
from rprobe import literal_cases, run

OUT = Path(__file__).resolve().parent / "probe_translate.jsonl"
K = 1000
OPS = {"nearest", "precede", "follow", "distance", "findOverlaps", "countOverlaps", "overlapsAny"}


def _shift_vec(v: dict) -> dict:
    v = dict(v)
    for k in ("int", "dbl"):
        if k in v:
            v[k] = [None if x is None else x + K for x in v[k]]
    return v


def translate(case: dict) -> dict:
    c = copy.deepcopy(case)
    for lit in c["inputs"].values():
        body = lit.get("IRanges") or lit.get("GRanges")
        if body is None:
            continue
        for k in ("start", "end"):
            if k in body:
                body[k] = _shift_vec(body[k])
    c["id"] += "#shifted"
    return c


def _circular_with_length(projected_inputs: dict) -> bool:
    """A range on a circular sequence of known length: translation is then not a symmetry."""
    used = {s for p in projected_inputs.values() if p.get("kind") == "GRanges" for s in p["seqnames"]}
    for p in projected_inputs.values():
        si = p.get("seqinfo") if p.get("kind") == "GRanges" else p if p.get("kind") == "Seqinfo" else None
        if si and any(c and n is not None and s in used
                      for s, c, n in zip(si["seqnames"], si["isCircular"], si["seqlengths"])):
            return True
    return False


def main() -> None:
    expected = {r["id"]: r for r in map(json.loads, (EVAL / "build" / "expected.jsonl").open())}
    resolved = {r["id"]: r for r in map(json.loads, (EVAL / "build" / "resolved.jsonl").open())}
    todo = []
    for cid, c in literal_cases().items():
        e = expected.get(cid)
        if e is None or c["op"] not in OPS or _circular_with_length(resolved[cid]["inputs"]):
            continue
        todo.append(translate(c))
    got = run(todo)
    n = 0
    with OUT.open("w") as fh:
        for tc in todo:
            cid = tc["id"].removesuffix("#shifted")
            r1 = got[tc["id"]].get("result")
            dep = None if r1 is None else r1 != expected[cid]["result"]  # None: shifted input not buildable
            n += bool(dep)
            fh.write(json.dumps({"id": cid, "translation_dependent": dep, "r_shifted": r1}) + "\n")
    print(f"{len(todo)} nearest/overlap cases re-run shifted by {K}; R's result changes in {n}")


if __name__ == "__main__":
    main()
