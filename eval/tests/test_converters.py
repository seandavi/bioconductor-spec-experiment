"""Converter rules (spec/CONVERTERS.md) on synthetic specs, plus reconciliation of the committed pool."""

import json

import pytest
from jsonschema import Draft202012Validator

from bioc_eval.common import CASES, SPEC, read_jsonl
from bioc_eval.convert import e1, e2, e4
from bioc_eval.convert.common import Builder, Unconvertible, as_integer, scalar_arg, untyped_vector, vec

CASE_SCHEMA = Draft202012Validator(json.loads((SPEC / "case.schema.json").read_text()))


def valid(case):
    assert [e.message for e in CASE_SCHEMA.iter_errors(case)] == []
    return case


# ---- shared literal rules
def test_untyped_vector_rules():
    assert untyped_vector([True, None], empty="int") == vec("lgl", [True, None])
    assert untyped_vector([1, 2], empty="int") == vec("int", [1, 2])
    assert untyped_vector([1, 2.5], empty="int") == vec("dbl", [1, 2.5])
    assert untyped_vector([1, 2 ** 31], empty="int") == vec("dbl", [1.0, float(2 ** 31)])
    assert untyped_vector(["a", None], empty="int") == vec("chr", ["a", None])
    assert untyped_vector([None, None], empty="int") == vec("lgl", [None, None])
    assert untyped_vector([], empty="int") == vec("int", [])
    assert untyped_vector([], empty="lgl") == vec("lgl", [])
    with pytest.raises(Unconvertible):
        untyped_vector([1, "a"], empty="int")
    with pytest.raises(Unconvertible):
        untyped_vector([], empty=None)


def test_scalar_arg():
    assert scalar_arg(3) == 3 and scalar_arg(True) is True and scalar_arg("s") == "s" and scalar_arg(2.5) == 2.5
    assert scalar_arg(2.0) == {"dbl": [2.0]}
    assert scalar_arg(2 ** 31) == {"dbl": [float(2 ** 31)]}
    assert scalar_arg(float("inf")) == {"dbl": ["Inf"]}


def test_as_integer_mirrors_r_for_exact_cases():
    assert as_integer([1, True, 2.9, None, -2.9]) == [1, 1, 2, None, -2]
    with pytest.raises(Unconvertible):
        as_integer(["3"])
    with pytest.raises(Unconvertible):
        as_integer([2 ** 31])


# ---- E1
def e1_case(spec, stem="x"):
    return valid(e1.convert_case(spec, stem, f"tests/differential/cases/{stem}.yaml"))


def test_e1_basic_call_and_refs():
    c = e1_case({"id": "a", "tags": ["findOverlaps"],
                 "inputs": {"q": {"IRanges": {"start": [1, 5], "width": [3, 0], "names": ["a", "b"], "mcols": {"m": [1, 2]}}}},
                 "call": {"fn": "findOverlaps", "args": ["$q", "$q"], "kwargs": {"type": "within", "maxgap": 2}}})
    assert c["id"] == "e1/x/a" and c["source"] == "e1" and c["op"] == "findOverlaps"
    assert c["inputs"]["q"]["IRanges"]["start"] == {"int": [1, 5]}
    assert c["inputs"]["q"]["IRanges"]["mcols"] == [["m", {"int": [1, 2]}]]
    assert c["call"]["kwargs"] == {"type": "within", "maxgap": 2}


def test_e1_redteam_source_by_file_or_tag():
    base = {"id": "a", "inputs": {"x": {"IRanges": {"start": [1], "width": [1]}}}, "call": {"fn": "start", "args": ["$x"]}}
    assert e1_case({**base, "tags": []}, "redteam_edges")["source"] == "e1-redteam"
    assert e1_case({**base, "tags": ["redteam"]}, "other")["source"] == "e1-redteam"
    assert e1_case({**base, "tags": []}, "other")["source"] == "e1"


def test_e1_null_is_lgl_na_and_NULL_literal_is_r_null_and_escape():
    c = e1_case({"id": "a", "tags": [], "inputs": {"x": {"IRanges": {"start": [1], "width": [1]}}},
                 "call": {"fn": "names<-", "args": ["$x"], "kwargs": {"value": {"NULL": []}}}})
    assert c["call"]["kwargs"]["value"] is None
    c = e1_case({"id": "b", "tags": [], "inputs": {"x": {"IRanges": {"start": [1], "width": [1]}}},
                 "call": {"fn": "shift", "args": ["$x", None]}})
    assert c["call"]["args"][1] == {"lgl": [None]}
    c = e1_case({"id": "c", "tags": [], "inputs": {"x": {"IRanges": {"start": [1], "width": [1]}}},
                 "call": {"fn": "resize", "args": ["$x", 2], "kwargs": {"fix": "$$start"}}})
    assert c["call"]["kwargs"]["fix"] == "$start"


def test_e1_identity_becomes_constructor_tree_with_hoisted_mcols():
    c = e1_case({"id": "i", "tags": [], "inputs": {"x": {"GRanges": {
        "seqnames": ["a"], "start": [1], "width": [2], "strand": ["+"],
        "seqinfo": {"seqnames": ["a", "b"], "seqlengths": [10, 20], "is_circular": [False, None]},
        "mcols": {"s": [1.5]}}}}, "call": {"fn": "identity", "args": ["$x"]}})
    assert c["op"] == "mcols<-"
    inner = c["call"]["args"][0]["$call"]
    assert inner["op"] == "GRanges"
    assert inner["kwargs"]["ranges"]["$call"]["op"] == "IRanges"
    assert inner["kwargs"]["seqinfo"]["$call"]["kwargs"]["isCircular"] == {"lgl": [False, None]}
    assert inner["kwargs"]["seqinfo"]["$call"]["kwargs"]["seqlengths"] == {"int": [10, 20]}
    df = c["call"]["kwargs"]["value"]["$ref"]
    assert c["inputs"][df]["DataFrame"]["columns"] == [["s", {"dbl": [1.5]}]]
    assert "x" not in c["inputs"]  # only referenced inputs are kept


def test_e1_object_literal_in_kwargs_is_hoisted():
    c = e1_case({"id": "h", "tags": [], "inputs": {"x": {"IRanges": {"start": [1], "width": [1]}}},
                 "call": {"fn": "c", "args": ["$x", {"IRanges": {"start": [4], "width": [1]}}]}})
    ref = c["call"]["args"][1]["$ref"]
    assert ref.startswith("_a") and "IRanges" in c["inputs"][ref]


def test_e1_raw_cases_and_unconvertibles():
    c = e1_case({"id": "rt-meta-gr-recycle-range-names", "tags": ["redteam"], "raw": {"r": "...", "py": "..."}})
    assert c["op"] == "GRanges" and c["call"]["args"][1]["$call"]["kwargs"] == {"width": 2, "names": "only"}
    with pytest.raises(Unconvertible) as ei:
        e1.convert_case({"id": "z", "raw": {"r": "1+1"}}, "x", "f")
    assert ei.value.reason == "raw-r"
    with pytest.raises(Unconvertible) as ei:
        e1.convert_case({"id": "z", "tags": [], "inputs": {"x": {"IRanges": {"start": [1]}, "width": [1]}},
                         "call": {"fn": "start", "args": ["$x"]}}, "x", "f")
    assert ei.value.reason == "literal:source-multi-key"
    with pytest.raises(Unconvertible) as ei:
        e1.convert_case({"id": "z", "tags": [], "inputs": {"x": [1, 2]}, "call": {"fn": "runValue", "args": ["$x"]}}, "x", "f")
    assert ei.value.reason == "op-out-of-vocab:runValue"


def test_e1_typed_empty_and_default_empty_int():
    c = e1_case({"id": "e", "tags": [], "inputs": {"x": {"IRanges": {"start": [], "width": {"int": []}}}},
                 "call": {"fn": "start", "args": ["$x"]}})
    assert c["inputs"]["x"]["IRanges"]["start"] == {"int": []}


# ---- E2
def test_e2_granges_defaults_made_explicit_and_fixture_closure():
    fx = {"si": {"$type": "Seqinfo", "seqnames": ["a", "b"], "seqlengths": {"$int": [5, None]}, "genome": "hg"},
          "unused": {"$type": "IRanges", "start": [1], "width": [1]},
          "g": {"$type": "GRanges", "seqnames": ["b", "a", "b"], "start": [1, 2, 3], "width": [1, 1, 1], "seqinfo": {"$ref": "si"}},
          "h": {"$type": "GRanges", "seqnames": ["z", "y", "z"], "start": [1, 2, 3], "width": [1, 1, 1]}}
    dec = e2.Decoder(fx)
    g = dec.fixture("g")["GRanges"]
    assert g["strand"] == vec("chr", ["*"] * 3) and g["seqinfo"] == {"$ref": "si"}
    h = dec.fixture("h")["GRanges"]
    assert h["seqinfo"] == {"seqnames": vec("chr", ["z", "y"])}  # Seqinfo(unique(seqnames)): appearance order
    assert dec.fixture("si")["Seqinfo"]["genome"] == vec("chr", ["hg"])


def test_e2_shorthand_and_tags():
    d = e2.Decoder({})
    assert d.shorthand([1, 2]) == vec("int", [1, 2]) and d.shorthand([1, 2.5]) == vec("dbl", [1, 2.5])
    assert d.shorthand([None, None]) == vec("lgl", [None, None])
    assert d.shorthand(["a"]) == vec("chr", ["a"])
    with pytest.raises(Unconvertible):
        d.shorthand([])
    with pytest.raises(Unconvertible):
        d.shorthand([1, "a"])
    assert d.vector({"$dbl": [1.0, "NaN", None], "names": ["a", "b", "c"]}) == vec("dbl", [1.0, "NaN", None], ["a", "b", "c"])
    with pytest.raises(Unconvertible) as ei:
        d.value_literal({"$type": "Hits", "from": {"$int": []}})
    assert ei.value.reason == "input-kind:Hits"


def test_e2_dframe_keeps_nrow_and_factor():
    d = e2.Decoder({})
    df = d.value_literal({"$type": "DFrame", "columns": {"f": {"$type": "Factor", "levels": {"$chr": ["a", "b"]}, "codes": {"$int": [1, None, 2]}}}})
    assert df["DataFrame"]["nrow"] == 3
    assert df["DataFrame"]["columns"][0][1] == {"Factor": {"levels": ["a", "b"], "codes": [1, None, 2]}}


# ---- E4
def test_e4_integer_coercion_names_and_refs():
    spec = {"id": "x", "clause": "C", "inputs": {
        "si": {"Seqinfo": {"seqnames": ["a", "b"], "seqlengths": [10.9, 20]}},
        "g": {"GRanges": {"seqnames": ["a"], "start": [1.7], "width": [True], "seqinfo": "si", "mcols": {"m": [1, None]}}}},
        "call": {"op": "nearest", "args": ["g", "g"], "kwargs": {"select": "all", "maxx": {"ref": "g"}, "v": [1, 2], "n": None}}}
    c = valid(e4.convert_case(spec, "granges", "overlaps", "cases/granges/overlaps.yaml"))
    assert c["id"] == "e4/granges-overlaps/x"
    g = c["inputs"]["g"]["GRanges"]
    assert g["start"] == {"int": [1]} and g["width"] == {"int": [1]} and g["seqinfo"] == {"$ref": "si"}
    assert c["inputs"]["si"]["Seqinfo"]["seqlengths"] == {"int": [10, 20]}
    assert c["call"]["kwargs"] == {"select": "all", "maxx": {"$ref": "g"}, "v": {"int": [1, 2]}, "n": None}


def test_e4_empty_untyped_list_is_logical():
    spec = {"id": "x", "clause": "C", "inputs": {"x": {"IRanges": {"start": [1], "width": [1]}}},
            "call": {"op": "shift", "args": ["x"], "kwargs": {"shift": []}}}
    assert e4.convert_case(spec, "iranges", "intra", "f")["call"]["kwargs"]["shift"] == {"lgl": []}


def test_e4_identity_non_object_is_out_of_vocab():
    spec = {"id": "x", "clause": "C", "inputs": {"v": {"int": [1, 2]}}, "call": {"op": "identity", "args": ["v"]}}
    with pytest.raises(Unconvertible) as ei:
        e4.convert_case(spec, "iranges", "edge", "f")
    assert ei.value.reason == "op-out-of-vocab:identity"


# ---- committed pool
def test_committed_pool_reconciles_and_validates():
    summary = json.loads((CASES / "conversion-summary.json").read_text())
    for suite, s in summary.items():
        assert s["read"] == s["converted"] + s["unconvertible"], suite
        assert s["converted"] == len(read_jsonl(CASES / f"{suite}.jsonl")), suite
        assert sum(s["unconvertible_by_reason"].values()) == s["unconvertible"]
    unconv = read_jsonl(CASES / "unconvertible.jsonl")
    assert len(unconv) == sum(s["unconvertible"] for s in summary.values())
    ids = set()
    for suite in ("e1", "e2", "e4", "heldout"):
        for case in read_jsonl(CASES / f"{suite}.jsonl"):
            valid(case)
            assert case["id"] not in ids
            ids.add(case["id"])
