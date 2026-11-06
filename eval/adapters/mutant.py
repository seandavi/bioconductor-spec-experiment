"""Self-test adapter: r_echo plus a deterministic perturbation on cases with int(hash[:8], 16) % 5 == 0.

Writes build/results/mutant.predicted.jsonl: the outcome the comparator must assign to every case.
Perturbation = first applicable entry of the 9-entry list, rotated by int(hash[8:16], 16) % 9.
Bump / drop are not applied to R3 cases (select="arbitrary" accepts other values), see NOTES.md.
"""

import copy
import json
from pathlib import Path

from adapter_base import AdapterBase, Unsupported

import fakeport

PERTURBATIONS = ["bump", "drop", "reverse_to", "swap_universe", "return_null",
                 "raise_call", "raise_build", "adapter_typeerror", "unsupported"]
DROP_KEYS = ("values", "start", "from", "lengths", "seqnames", "elements", "columns", "mcols", "rownames", "to", "width")


def first_leaf(node, path=()):
    """Path to the first scalar leaf (None excluded), skipping `kind` / `type` keys."""
    if isinstance(node, dict):
        for k, v in node.items():
            if k in ("kind", "type"):
                continue
            if (r := first_leaf(v, path + (k,))) is not None:
                return r
    elif isinstance(node, list):
        for i, v in enumerate(node):
            if (r := first_leaf(v, path + (i,))) is not None:
                return r
    elif node is not None:
        return path
    return None


def bump(value):
    path = first_leaf(value)
    node = value
    for k in path[:-1]:
        node = node[k]
    v = node[path[-1]]
    if isinstance(v, bool):
        new = not v
    elif isinstance(v, int):
        new = v + 1
    elif isinstance(v, float):
        new = v + 1.0 if v + 1.0 != v else v * 2.0
    else:
        new = v + "x"
    node[path[-1]] = new
    return value


def drop_key(value):
    for k in DROP_KEYS:
        if isinstance(value.get(k), list) and value[k]:
            return k
    return None


def multi_hit_query(hits):
    """(start, stop) of the first run of equal `from` with >= 2 hits, else None."""
    f = hits["from"]
    i = 0
    while i < len(f):
        j = i
        while j < len(f) and f[j] == f[i]:
            j += 1
        if j - i >= 2:
            return i, j
        i = j
    return None


def universe_hits(exp):
    """Per query: distinct allowed picks (R's select="all" hits plus R's own pick), sorted."""
    by: dict[int, set[int]] = {}
    for f, t in zip(exp["universe"]["from"], exp["universe"]["to"]):
        by.setdefault(f, set()).add(t)
    for i, v in enumerate(exp["result"]["values"], start=1):
        if v is not None:
            by.setdefault(i, set()).add(v)
    return {q: sorted(v) for q, v in by.items()}


class Adapter(AdapterBase):
    port = "mutant"
    port_packages = ("fakeport",)

    def setup(self, ctx: dict) -> None:
        self.expected = {}
        for line in Path(ctx["expected_path"]).read_text().splitlines():
            rec = json.loads(line)
            if rec["status"] == "ok":
                self.expected[rec["id"]] = rec
        self.predicted: dict[str, tuple[str, str | None]] = {}
        self.out_path = Path(ctx["out_path"]).with_name("mutant.predicted.jsonl")

    # ---- planning -------------------------------------------------------------------
    def applicable(self, exp: dict, name: str) -> bool:
        res = exp["result"]
        r_error = res["kind"] == "error"
        r3 = "universe" in exp and res["kind"] == "vector"
        if name == "bump":
            return not r_error and not r3 and first_leaf(res) is not None
        if name == "drop":
            return not r_error and not r3 and drop_key(res) is not None
        if name == "reverse_to":
            return res["kind"] == "Hits" and bool(exp.get("sorted_by_query")) and multi_hit_query(res) is not None
        if name == "swap_universe":
            return r3 and any(len(v) >= 2 for v in universe_hits(exp).values())
        if name == "return_null":
            return r_error
        if name == "raise_call":
            return not r_error
        if name == "raise_build":
            return r_error
        return True  # adapter_typeerror, unsupported

    def plan(self, case: dict) -> tuple[str | None, str]:
        h = case["hash"]
        if int(h[:8], 16) % 5 != 0:
            return None, "pass"
        exp = self.expected[case["id"]]
        start = int(h[8:16], 16) % len(PERTURBATIONS)
        for k in range(len(PERTURBATIONS)):
            name = PERTURBATIONS[(start + k) % len(PERTURBATIONS)]
            if self.applicable(exp, name):
                predicted = {"bump": "mismatch", "drop": "mismatch", "reverse_to": "pass", "swap_universe": "pass",
                             "return_null": "mismatch", "raise_call": "port_error_unexpected",
                             "raise_build": "port_error_unexpected", "adapter_typeerror": "adapter_error",
                             "unsupported": "unsupported"}[name]
                return name, predicted
        raise AssertionError("unreachable: adapter_typeerror always applies")

    # ---- execution --------------------------------------------------------------------
    def build_phase(self, exp: dict, name: str | None) -> None:
        if name == "raise_build":
            fakeport.boom("mutant: raise during build")
        if name == "adapter_typeerror":
            raise TypeError("mutant: translation error in adapter code")
        if name == "unsupported":
            raise Unsupported("mutant: declared unsupported")

    def call_phase(self, exp: dict, name: str | None) -> dict:
        res = copy.deepcopy(exp["result"])
        if res["kind"] == "error":
            if name == "return_null":
                return {"kind": "NULL"}
            fakeport.boom("R error: " + res["message"][:100])
        if name == "raise_call":
            fakeport.boom("mutant: raise during call")
        if name == "bump":
            return bump(res)
        if name == "drop":
            res[drop_key(res)].pop()
            return res
        if name == "reverse_to":
            i, j = multi_hit_query(res)
            res["to"][i:j] = res["to"][i:j][::-1]
            return res
        if name == "swap_universe":
            for q, hits in sorted(universe_hits(exp).items()):
                if len(hits) >= 2:
                    cur = res["values"][q - 1]
                    res["values"][q - 1] = next(t for t in hits if t != cur)
                    return res
        return res

    def run_case(self, case: dict, stage) -> dict:
        exp = self.expected[case["id"]]
        name, predicted = self.plan(case)
        self.predicted[case["id"]] = (predicted, name)
        stage("build", self.build_phase, exp, name)
        return {"result": stage("call", self.call_phase, exp, name)}

    def finish(self) -> None:
        with self.out_path.open("w") as fh:
            for cid in sorted(self.predicted):
                predicted, name = self.predicted[cid]
                fh.write(json.dumps({"id": cid, "predicted": predicted, "perturbation": name}) + "\n")
