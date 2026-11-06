"""Oracle: ten hand-checked cases (answers derived by hand from the R documentation, not from R output),
a determinism check, and the projection of R errors / exclusions."""

import shutil

import pytest

from bioc_eval.common import load_ops
from bioc_eval.oracle import run_r, universe_call

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="needs Rscript")


def ir(start, width, **kw):
    return {"IRanges": {"start": {"int": start}, "width": {"int": width}, **kw}}


def call(op, *args, **kwargs):
    return {"op": op, "args": list(args), "kwargs": kwargs}


def ref(n):
    return {"$ref": n}


def vecp(vals, type_="int", names=None):
    return {"kind": "vector", "type": type_, "values": vals, "names": names}


def irp(start, width):
    return {"kind": "IRanges", "start": start, "width": width, "names": None, "mcols": []}


HAND = [
    # 1 reduce merges (1..3)+(2..5), keeps (10..12)
    ("reduce", {"x": ir([1, 2, 10], [3, 4, 3])}, call("reduce", ref("x")), irp([1, 10], [5, 3])),
    # 2 gaps between 2..4 and 10..12 within 1..15
    ("gaps", {"x": ir([2, 10], [3, 3])}, call("gaps", ref("x"), start=1, end=15), irp([1, 5, 13], [1, 5, 3])),
    # 3 3..4 lies within 1..10
    ("findOverlaps", {"q": ir([3], [2]), "s": ir([1], [10])}, call("findOverlaps", ref("q"), ref("s"), type="within"),
     {"kind": "Hits", "from": [1], "to": [1], "nLnode": 1, "nRnode": 1}),
    # 4 shift by 3
    ("shift", {"x": ir([1, 5], [2, 2])}, call("shift", ref("x"), 3), irp([4, 8], [2, 2])),
    # 5 5..7 resized to width 1 fixing the end -> 7..7
    ("resize", {"x": ir([5], [3])}, call("resize", ref("x"), 1, fix="end"), irp([7], [1])),
    # 6 coverage of 1..3 and 3..5
    ("coverage", {"x": ir([1, 3], [3, 3])}, call("coverage", ref("x")),
     {"kind": "Rle", "type": "int", "values": [1, 2, 1], "lengths": [2, 1, 2]}),
    # 7 minus-strand flank(width 5) is upstream = past the end: 20..24
    ("flank", {"g": {"GRanges": {"seqnames": {"chr": ["chr1"]}, "start": {"int": [10]}, "width": {"int": [10]}, "strand": {"chr": ["-"]}}}},
     call("flank", ref("g"), 5), None),
    # 8 negative width is an error
    ("IRanges", {}, call("IRanges", start={"int": [1]}, width={"int": [-1]}), {"kind": "error"}),
    # 9 order by start then width
    ("order", {"x": ir([3, 1, 2], [1, 1, 1])}, call("order", ref("x")), vecp([2, 3, 1])),
    # 10 countOverlaps: q1 1..5 hits s1 3..4 and s2 4..5; q2 10..11 hits nothing
    ("countOverlaps", {"q": ir([1, 10], [5, 2]), "s": ir([3, 4, 20], [2, 2, 2])}, call("countOverlaps", ref("q"), ref("s")),
     vecp([2, 0])),
]


def case(i, op, inputs, c):
    return {"id": f"t/{i}", "source": "heldout", "op": op, "inputs": inputs, "call": c}


def test_hand_checked_cases():
    ops = load_ops()
    cases = [case(i, op, inp, c) for i, (op, inp, c, _) in enumerate(HAND)]
    got = run_r(cases, ops, workers=1)
    for i, (op, _, _, want) in enumerate(HAND):
        rec = got[f"t/{i}"]
        assert rec["status"] == "ok", rec
        res = rec["result"]
        if op == "flank":
            assert (res["kind"], res["start"], res["width"], res["strand"]) == ("GRanges", [20], [5], ["-"])
        elif want == {"kind": "error"}:
            assert res["kind"] == "error" and "non-negative" in res["message"]
        else:
            assert res == want, (op, res)


def test_r_input_error_excludes_and_unprojectable_excludes():
    ops = load_ops()
    bad = case(0, "start", {"x": {"IRanges": {"start": {"int": [1]}, "width": {"int": [-5]}}}}, call("start", ref("x")))
    unproj = case(1, "split", {}, call("list", {"int": [1]}))  # `list` is not in the vocabulary but R can run it
    got = run_r([bad, unproj], ops, workers=1)
    assert got["t/0"]["status"] == "excluded" and got["t/0"]["reason"] == "r_input_error"
    assert got["t/1"]["status"] == "excluded" and got["t/1"]["reason"].startswith("r_unprojectable:")


def test_hits_class_and_universe_recorded_for_arbitrary():
    ops = load_ops()
    inputs = {"q": ir([1, 5], [3, 3]), "s": ir([2, 6, 6], [10, 2, 0])}
    c = case(0, "findOverlaps", inputs, call("findOverlaps", ref("q"), ref("s"), select="arbitrary"))
    u = universe_call(c, ops)
    assert u["kwargs"]["select"] == "all" and u["args"] == []
    h = case(1, "findOverlaps", inputs, call("findOverlaps", ref("q"), ref("s")))
    assert universe_call(h, ops) is None
    n = case(2, "nearest", inputs, call("nearest", ref("q"), ref("s")))
    assert universe_call(n, ops)["kwargs"]["select"] == "all"
    got = run_r([c, h], ops, workers=1)
    assert got["t/0"]["universe"]["kind"] == "Hits"
    assert got["t/1"]["sorted_by_query"] is True and "universe" not in got["t/1"]


def test_select_given_positionally_is_found():
    ops = load_ops()
    inputs = {"q": ir([1], [3]), "s": ir([2], [3])}
    c = case(0, "findOverlaps", inputs, call("findOverlaps", ref("q"), ref("s"), -1, 0, "any", "arbitrary"))
    assert universe_call(c, ops)["kwargs"]["select"] == "all"


def test_second_run_identical():
    ops = load_ops()
    cases = [case(i, op, inp, c) for i, (op, inp, c, _) in enumerate(HAND)]
    assert run_r(cases, ops, workers=1) == run_r(cases, ops, workers=2)
