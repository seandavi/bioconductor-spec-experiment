"""Random differential testing: R vs ranges_e2 on generated cases (nothing is committed).

    uv run python -m ranges_e2.conformance.fuzz --seed 1 --n 2000 [--min-start -2] [--show 5]

Writes a throw-away `fuzz` case group into a temp copy of conformance/, has R generate goldens,
runs the Python runner on it and prints the divergences (minus the documented empty-subject quirk
of GRanges precede/follow). A divergence is either a Python bug,
an R quirk worth a named clause + case in the real suite, or out-of-domain input (see README).
"""

from __future__ import annotations

import argparse
import json
import random
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

from .codec import decode, encode
from .ops import run_op
from .runner import CONFORMANCE_DIR, load_cases, load_expected, status_of

SEQINFO = {
    "$type": "Seqinfo",
    "seqnames": ["c1", "c2", "c3"],
    "seqlengths": {"$int": [40, 30, None]},
    "isCircular": [False, False, False],
    "genome": "g",
}
OPS = (
    "shift narrow resize flank promoters range reduce gaps disjoin union intersect setdiff punion pintersect "
    "psetdiff find count any subset nearest precede follow distance coverage order sort rank dup unique trim"
).split()


def _tag(values: list[Any], kind: str = "$int") -> dict[str, Any]:
    return {kind: values}


class Generator:
    def __init__(self, seed: int, min_start: int) -> None:
        self.rng = random.Random(seed)
        self.min_start = min_start
        self.cases: list[dict[str, Any]] = []

    def ranges(self, n: int | None = None, genomic: bool = False) -> dict[str, Any]:
        r = self.rng
        n = r.randint(0, 7) if n is None else n
        obj: dict[str, Any] = {
            "$type": "GRanges" if genomic else "IRanges",
            "start": _tag([r.randint(self.min_start, 25) for _ in range(n)]),
            "width": _tag([r.choice([0, 0, 1, 2, 3, 5, 8]) for _ in range(n)]),
        }
        if genomic:
            obj["seqnames"] = _tag([r.choice(["c1", "c2", "c3"]) for _ in range(n)], "$chr")
            obj["strand"] = _tag([r.choice(["+", "-", "*"]) for _ in range(n)], "$chr")
            obj["seqinfo"] = SEQINFO
        return obj

    def add(self, op: str, args: list[Any], kwargs: dict[str, Any], tag: str) -> None:
        self.cases.append({"id": f"fz.{len(self.cases):05d}.{tag}.{op}", "clause": "FZ", "op": op, "args": args, "kwargs": kwargs})

    def one(self) -> None:
        r = self.rng
        genomic = r.random() < 0.5
        x, y = self.ranges(genomic=genomic), self.ranges(genomic=genomic)
        n = len(x["start"]["$int"])
        tag = "gr" if genomic else "ir"
        ign = {"ignore.strand": r.choice([True, False])} if genomic and r.random() < 0.4 else {}
        pick = r.choice(OPS)

        def i(lo: int, hi: int) -> dict[str, Any]:
            return _tag([r.randint(lo, hi)])

        def add(op: str, args: list[Any], kw: dict[str, Any]) -> None:
            self.add(op, args, kw, tag)

        if pick == "shift":
            add("shift", [x, i(-5, 5)], {})
        elif pick == "narrow":
            add("narrow", [x], r.choice([{"start": i(-4, 4)}, {"end": i(-4, 4)}, {"start": i(1, 3), "width": i(0, 3)}, {"start": i(-5, 3), "end": i(-3, 5)}]))
        elif pick == "resize":
            add("resize", [x, i(0, 9)], {"fix": r.choice(["start", "end", "center"]), **ign})
        elif pick == "flank":
            add("flank", [x, i(-4, 6)], {"start": r.choice([True, False]), "both": r.choice([True, False]), **ign})
        elif pick == "promoters":
            add("promoters", [x], {"upstream": i(0, 8), "downstream": i(0, 8)})
        elif pick == "range":
            add("range", [x], dict(ign))
        elif pick == "reduce":
            add("reduce", [x], {"min.gapwidth": r.randint(0, 3), "drop.empty.ranges": r.choice([True, False]), "with.revmap": r.choice([True, False]), **ign})
        elif pick == "gaps":
            kw: dict[str, Any] = dict(ign)
            if r.random() < 0.5:
                kw["start"] = i(-3, 10)
            if r.random() < 0.5:
                kw["end"] = i(0, 35)
            add("gaps", [x], kw)
        elif pick == "disjoin":
            add("disjoin", [x], {"with.revmap": r.choice([True, False]), **ign})
        elif pick in ("union", "intersect", "setdiff"):
            add(pick, [x, y], dict(ign))
        elif pick in ("punion", "pintersect", "psetdiff"):
            kw = dict(ign)
            if pick == "punion":
                kw["fill.gap"] = r.choice([True, False])
            if pick == "pintersect" and not genomic:
                kw["resolve.empty"] = r.choice(["none", "max.start", "start.x"])
            add(pick, [x, self.ranges(n, genomic)], kw)
        elif pick in ("find", "count", "any", "subset"):
            t = r.choice(["any", "start", "end", "within", "equal"])
            kw = {"type": t, **ign}
            if t == "any":
                if r.random() < 0.5:
                    kw["maxgap"] = r.randint(-1, 4)
                elif r.random() < 0.5:
                    kw["minoverlap"] = r.randint(0, 3)
            else:
                kw["maxgap"], kw["minoverlap"] = r.randint(-1, 4), r.randint(0, 2)
            if pick == "find":
                add("findOverlaps", [x, y], {**kw, "select": r.choice(["all", "first", "last"])})
            elif pick == "count":
                add("countOverlaps", [x, y], kw)
            elif pick == "any":
                add("overlapsAny", [x, y], kw)
            else:
                add("subsetByOverlaps", [x, y], {**kw, "invert": r.choice([True, False])})
        elif pick == "nearest":
            add("nearest", [x, y], {"select": "all", **ign})  # 'arbitrary' is not a contract on ties
        elif pick == "precede":
            add("precede", [x, y], {"select": r.choice(["first", "all"]), **ign})
        elif pick == "follow":
            add("follow", [x, y], {"select": r.choice(["last", "all"]), **ign})
        elif pick == "distance":
            add("distance", [x, self.ranges(n if r.random() < 0.7 else 1, genomic)], dict(ign))
        elif pick == "coverage":
            kw = {}
            if r.random() < 0.4:
                kw["shift"] = i(0, 4)
            if r.random() < 0.4 and not genomic:
                kw["width"] = i(0, 40)
            if r.random() < 0.4:
                kw["weight"] = r.choice([i(1, 3), {"$dbl": [r.choice([0.5, 1.5, 2.0])]}])
            add("coverage", [x], kw)
        elif pick == "order":
            add("order", [x], {"decreasing": r.choice([True, False])})
        elif pick == "sort":
            add("sort", [x], {"decreasing": r.choice([True, False]), **ign})
        elif pick == "rank":
            add("rank", [x], {"ties.method": r.choice(["average", "first", "last", "min", "max"]), **ign})
        elif pick in ("dup", "unique"):
            add("duplicated" if pick == "dup" else "unique", [x], {"fromLast": r.choice([True, False])})
        elif pick == "trim" and genomic:
            add("trim", [x], {})


def _known_quirk(case: Any) -> bool:
    """GRanges precede/follow with an empty query or subject return Hits in R (see xfail.json)."""
    if case.op not in ("precede", "follow") or case.args[0].get("$type") != "GRanges":
        return False
    return any(len(a["start"]["$int"]) == 0 for a in case.args[:2])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--min-start", type=int, default=1, help="smallest range start; <1 leaves R's supported domain")
    ap.add_argument("--show", type=int, default=5, help="print this many divergences in detail")
    args = ap.parse_args()
    gen = Generator(args.seed, args.min_start)
    for _ in range(args.n):
        gen.one()
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp) / "conformance"
        (base / "cases").mkdir(parents=True)
        (base / "expected").mkdir()
        (base / "reports").mkdir()
        shutil.copytree(CONFORMANCE_DIR / "r", base / "r")
        doc = {"group": "fuzz", "compare": "hits-query-grouped", "clauses": {"FZ": "random"}, "fixtures": {}, "cases": gen.cases}
        (base / "cases" / "fuzz.json").write_text(json.dumps(doc))
        subprocess.run(["Rscript", str(base / "r" / "run.R"), "generate", "fuzz"], check=True, capture_output=True)
        expected = load_expected(base)
        bad = [(c, *status_of(c, expected, {})) for c in load_cases(base)]
        bad = [b for b in bad if b[1] != "pass" and not _known_quirk(b[0])]
        print(f"{len(gen.cases)} cases, {len(bad)} divergences", dict(Counter(c.op for c, _, _ in bad)))
        for c, status, detail in bad[: args.show]:
            try:
                got: Any = encode(run_op(c.op, [decode(a) for a in c.args], {k: decode(v) for k, v in c.kwargs.items()}))
            except Exception as e:  # noqa: BLE001
                got = f"{type(e).__name__}: {e}"
            exp = expected[c.id]
            print(f"\n{c.id} [{status}] {detail}\n  args={json.dumps(c.args)[:300]}\n  kwargs={c.kwargs}")
            print(f"  python={json.dumps(got)[:300]}\n  R     ={json.dumps(exp.get('value', exp.get('error')))[:300]}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
