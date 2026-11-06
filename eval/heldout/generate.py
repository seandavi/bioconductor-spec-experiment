"""Seeded heldout case generator (stdlib only). Scope: TARGET.md. Never runs R.

    python3 heldout/generate.py            # writes cases/heldout.jsonl, heldout/COVERAGE.md, heldout/seeds.txt
    python3 heldout/generate.py --check    # regenerate in memory, fail if committed output differs

Per-case rng: seed = sha256("bioc-spec-eval/heldout/v1/<method>/<k>")[:8] (gen_core.seed_for).
"""

from __future__ import annotations

import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import gen_ops as O  # noqa: E402
from gen_core import G, seed_for  # noqa: E402

SPEC = HERE.parent / "spec"
OUT = HERE.parent / "cases" / "heldout.jsonl"
COVERAGE = HERE / "COVERAGE.md"
SEEDS = HERE / "seeds.txt"

BOTH = ("IRanges", "GRanges")
GR = ("GRanges",)
BASE_N = 30
BIG_N = 60
BIG = {"findOverlaps", "nearest", "reduce", "coverage"}

# method -> (classes, builder).  Method names are ops.json `method`.
REGISTRY: dict[str, tuple[tuple[str, ...], object]] = {
    "IRanges": (("IRanges",), O.op_IRanges),
    "GRanges": (GR, O.op_GRanges),
    "Seqinfo": (("Seqinfo",), O.op_Seqinfo),
    "start": (BOTH, O.op_plain_accessor), "end": (BOTH, O.op_plain_accessor),
    "width": (BOTH, O.op_plain_accessor), "names": (BOTH, O.op_plain_accessor),
    "length": (BOTH, O.op_plain_accessor),
    "seqnames": (GR, O.op_plain_accessor), "strand": (GR, O.op_plain_accessor),
    "ranges": (BOTH, O.op_ranges), "mcols": (BOTH, O.op_mcols),
    "set_names": (BOTH, O.op_set_names), "set_mcols": (BOTH, O.op_set_mcols),
    "seqinfo": (GR, O.op_seqinfo_getter),
    "seqlevels": (GR + ("Seqinfo",), O.op_seqinfo_getter),
    "seqlengths": (GR + ("Seqinfo",), O.op_seqinfo_getter),
    "isCircular": (GR + ("Seqinfo",), O.op_seqinfo_getter),
    "genome": (GR + ("Seqinfo",), O.op_seqinfo_getter),
    "set_seqinfo": (GR, O.op_set_seqinfo),
    "set_seqlevels": (GR, O.op_set_seqlevels),
    "set_seqlengths": (GR + ("Seqinfo",), O.op_set_seqlengths),
    "set_isCircular": (GR + ("Seqinfo",), O.op_set_isCircular),
    "set_genome": (GR + ("Seqinfo",), O.op_set_genome),
    "merge": (("Seqinfo",), O.op_merge),
    "trim": (GR, O.op_trim),
    "subset": (BOTH, O.op_subset), "c": (BOTH, O.op_c),
    "shift": (BOTH, O.op_shift), "narrow": (BOTH, O.op_narrow), "resize": (BOTH, O.op_resize),
    "flank": (BOTH, O.op_flank), "promoters": (BOTH, O.op_promoters),
    "range": (BOTH, O.op_range), "reduce": (BOTH, O.op_reduce), "gaps": (BOTH, O.op_gaps),
    "disjoin": (BOTH, O.op_disjoin),
    "union": (BOTH, O.op_setop), "intersect": (BOTH, O.op_setop), "setdiff": (BOTH, O.op_setop),
    "punion": (BOTH, O.op_parallel), "pintersect": (BOTH, O.op_parallel), "psetdiff": (BOTH, O.op_parallel),
    "findOverlaps": (BOTH, O.op_findOverlaps), "countOverlaps": (BOTH, O.op_countOverlaps),
    "overlapsAny": (BOTH, O.op_overlapsAny), "subsetByOverlaps": (BOTH, O.op_subsetByOverlaps),
    "nearest": (BOTH, O.op_nearest_like), "precede": (BOTH, O.op_nearest_like),
    "follow": (BOTH, O.op_nearest_like), "distance": (BOTH, O.op_distance),
    "coverage": (BOTH, O.op_coverage),
    "order": (BOTH, O.op_order), "sort": (BOTH, O.op_sort), "rank": (BOTH, O.op_rank),
    "duplicated": (BOTH, O.op_duplicated), "unique": (BOTH, O.op_unique),
}
# constructors need more cases: any two of start/end/width, error combos
COUNT_OVERRIDE = {"IRanges": 90, "GRanges": 90, "Seqinfo": 60}

FEATURES = ["empty", "zero_width", "names", "mcols", "strand", "seqinfo", "na_arg"]
EXTRA = ["seqinfo_unused", "seqinfo_unsorted", "seqinfo_lengths", "seqinfo_circular", "seqinfo_genome",
         "c_merge_seqinfo", "trim_seqlengths", "order_cross_seqlevel_strand", "revmap", "invalid", "extreme"]


# ------------------------------------------------------------------ feature detection

def _vec_has_na(node) -> bool:
    for key in ("int", "dbl", "chr", "lgl"):
        if isinstance(node, dict) and key in node:
            return any(v is None or v == "NaN" for v in node[key])
    return False


def _walk_args(call: dict):
    for a in [*call["args"], *call["kwargs"].values()]:
        yield a
        if isinstance(a, dict) and "$call" in a:
            yield from _walk_args(a["$call"])


def _literal_features(lit: dict) -> set[str]:
    f: set[str] = set()
    if "IRanges" in lit or "GRanges" in lit:
        gr = "GRanges" in lit
        body = lit["GRanges" if gr else "IRanges"]
        lens = [len(body[k][t]) for k in ("start", "end", "width") if k in body for t in body[k] if t in ("int", "dbl")]
        if gr:
            lens.append(len(body["seqnames"]["chr"]))
        if lens and min(lens) == 0:
            f.add("empty")
        if "width" in body and 0 in body["width"].get("int", []):
            f.add("zero_width")
        if "start" in body and "end" in body:
            if any(e - s == -1 for s, e in zip(body["start"]["int"], body["end"]["int"])):
                f.add("zero_width")
        if "start" in body and "width" not in body and "end" in body:
            pass
        if "end" in body and "width" in body:
            if 0 in body["width"]["int"]:
                f.add("zero_width")
        if "names" in body:
            f.add("names")
        if body.get("mcols"):
            f.add("mcols")
        if gr:
            st = body.get("strand", {}).get("chr", [])
            if len(set(st)) > 1:
                f.add("strand")
            si = body.get("seqinfo")
            if isinstance(si, dict) and "seqnames" in si:
                f |= _seqinfo_features(si)
                if set(body["seqnames"]["chr"]) < set(si["seqnames"]["chr"]):
                    f.add("seqinfo_unused")
    elif "Seqinfo" in lit:
        f |= _seqinfo_features(lit["Seqinfo"])
    elif "DataFrame" in lit:
        if lit["DataFrame"]["columns"]:
            f.add("mcols")
    return f


def _seqinfo_features(si: dict) -> set[str]:
    f: set[str] = set()
    lv = si["seqnames"]["chr"]
    if len(lv) > 1:
        f.add("seqinfo")
    if lv != sorted(lv):
        f.add("seqinfo_unsorted")
    if any(v is not None for v in si.get("seqlengths", {}).get("int", [])):
        f |= {"seqinfo", "seqinfo_lengths"}
    if any(si.get("isCircular", {}).get("lgl", [])):
        f |= {"seqinfo", "seqinfo_circular"}
    if si.get("genome", {}).get("chr"):
        f |= {"seqinfo", "seqinfo_genome"}
    return f


def case_features(inputs: dict, call: dict, tags: list[str]) -> list[str]:
    f = set(tags)
    for lit in inputs.values():
        f |= _literal_features(lit)
    args = list(_walk_args(call))
    if any(_vec_has_na(a) for a in args):
        f.add("na_arg")
    if "ignore.strand" in call["kwargs"]:
        f.add("strand")
    ordered = [t for t in FEATURES + EXTRA if t in f]
    return ordered


def reachable_inputs(inputs: dict, call: dict) -> dict:
    """Keep only inputs the call references (transitively), in original order."""
    need: set[str] = set()

    def refs(node):
        if isinstance(node, dict):
            if "$ref" in node and len(node) == 1:
                need.add(node["$ref"])
            else:
                for v in node.values():
                    refs(v)
        elif isinstance(node, list):
            for v in node:
                refs(v)

    refs(call)
    changed = True
    while changed:
        before = len(need)
        for nm in list(need):
            if nm in inputs:
                refs(inputs[nm])
        changed = len(need) != before
    missing = need - set(inputs)
    if missing:
        raise AssertionError(f"dangling $ref {sorted(missing)}")
    return {k: v for k, v in inputs.items() if k in need}


# ------------------------------------------------------------------ generation

def build_case(entries: dict, method: str, cls: str, k: int, builder) -> dict:
    seed = seed_for(method, k)
    rng = random.Random(seed)
    g = G(rng, entries[method], cls, k)
    call = builder(g)
    inputs = reachable_inputs(g.inputs, call)
    tags = case_features(inputs, call, g.tags)
    return {
        "id": f"heldout/{method}/{cls}/{k}",
        "source": "heldout",
        "origin": {"suite": "heldout", "file": "heldout/generate.py", "id": f"{method}/{cls}/{k}", "tags": tags},
        "op": entries[method]["op"],
        "inputs": inputs,
        "call": call,
    }


def generate() -> tuple[list[dict], list[int]]:
    ops = json.loads((SPEC / "ops.json").read_text())["ops"]
    entries = {o["method"]: o for o in ops}
    missing = sorted(set(entries) - set(REGISTRY))
    extra = sorted(set(REGISTRY) - set(entries))
    assert not missing and not extra, f"registry mismatch: missing {missing}, extra {extra}"
    cases, seeds = [], []
    for method in entries:  # ops.json order
        classes, builder = REGISTRY[method]
        per_class = COUNT_OVERRIDE.get(method, BIG_N if method in BIG else BASE_N)
        for ci, cls in enumerate(classes):
            for j in range(per_class):
                k = ci * per_class + j
                cases.append(build_case(entries, method, cls, k, builder))
                seeds.append(seed_for(method, k))
    cases.sort(key=lambda c: c["id"])
    return cases, seeds


def coverage_report(cases: list[dict], families: dict[str, str]) -> tuple[str, list[str]]:
    fam_feat: dict[str, Counter] = defaultdict(Counter)
    fam_n: Counter = Counter()
    opcls: Counter = Counter()
    for c in cases:
        fam = families[c["op"]]
        fam_n[fam] += 1
        _, method, cls, _ = c["id"].split("/")
        opcls[(c["op"], cls)] += 1
        for t in c["origin"]["tags"]:
            fam_feat[fam][t] += 1
    fams = list(dict.fromkeys(families.values()))
    problems = []
    lines = ["# Heldout coverage", "", f"{len(cases)} cases. Cells = cases carrying the feature (required: >= 5 per family for the first seven).", ""]
    cols = FEATURES + EXTRA
    lines.append("| family | n | " + " | ".join(cols) + " |")
    lines.append("|---|---|" + "---|" * len(cols))
    for fam in fams:
        row = [str(fam_feat[fam][t]) for t in cols]
        lines.append(f"| {fam} | {fam_n[fam]} | " + " | ".join(row) + " |")
        for t in FEATURES:
            if fam_feat[fam][t] < 5:
                problems.append(f"{fam}: feature {t} has {fam_feat[fam][t]} < 5")
    tot = Counter()
    for fam in fams:
        tot.update(fam_feat[fam])
    for t in ("c_merge_seqinfo", "trim_seqlengths", "order_cross_seqlevel_strand"):
        if tot[t] < 5:
            problems.append(f"{t} has {tot[t]} < 5")
    lines += ["", "Legend: empty = a range input of length 0; zero_width = a width-0 range; names; mcols; strand = mixed",
              "`+ - *` or `ignore.strand` given; seqinfo = multi-level / lengths / circular / genome; na_arg = NA in an",
              "atomic argument; invalid = a deliberately invalid argument value was injected; extreme = coordinates near int32 limits.",
              "", "## Cases per (op, class)", "", "| op | class | n |", "|---|---|---|"]
    for (op, cls), n in sorted(opcls.items()):
        lines.append(f"| {op} | {cls} | {n} |")
    inv = sum(1 for c in cases if "invalid" in c["origin"]["tags"])
    ext = sum(1 for c in cases if "extreme" in c["origin"]["tags"])
    lines += ["", f"Invalid-injected: {inv} ({100 * inv / len(cases):.1f} %); extreme-coordinate: {ext} ({100 * ext / len(cases):.1f} %)."]
    return "\n".join(lines) + "\n", problems


def dumps(cases: list[dict]) -> str:
    return "".join(json.dumps(c, ensure_ascii=False, separators=(",", ":")) + "\n" for c in cases)


def main() -> int:
    ops = json.loads((SPEC / "ops.json").read_text())["ops"]
    families = {o["op"]: o["family"] for o in ops}
    cases, seeds = generate()
    text = dumps(cases)
    report, problems = coverage_report(cases, families)
    if "--check" in sys.argv:
        ok = OUT.read_text() == text and COVERAGE.read_text() == report
        print("heldout output up to date" if ok else "heldout output DIFFERS from committed")
        return 0 if ok and not problems else 1
    OUT.write_text(text)
    COVERAGE.write_text(report)
    SEEDS.write_text("".join(f"{s}\n" for s in sorted(set(seeds))))
    print(report.split("## Cases per")[0])
    for p in problems:
        print("COVERAGE PROBLEM:", p)
    print(f"{len(cases)} cases -> {OUT}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
