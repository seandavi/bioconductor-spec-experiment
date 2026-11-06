"""Adapter helpers: RVector, call evaluation, p_* projection builders."""

import math

import numpy as np
import pytest

from adapter_base import (
    AdapterBase, RVector, Unsupported, UnknownArgument, decode_arg_literal, p_dataframe, p_factor, p_granges,
    p_hits, p_integerlist, p_iranges, p_null, p_rle, p_rlelist, p_seqinfo, p_vector, p_vector_from,
)
from binding import BindError


def test_decode_arg_literal():
    assert decode_arg_literal(None) is None
    assert decode_arg_literal(True).type == "lgl"
    assert decode_arg_literal(3).type == "int"
    assert decode_arg_literal(2.5).type == "dbl"
    assert decode_arg_literal("s").scalar() == "s"
    v = decode_arg_literal({"dbl": ["NaN", None, 1.0], "names": ["a", "b", "c"]})
    assert v.names == ["a", "b", "c"]
    lst = v.list()
    assert math.isnan(lst[0]) and lst[1] is None and lst[2] == 1.0


def test_rvector_numpy_views():
    v = RVector("int", [1, None, 3])
    with pytest.raises(ValueError):
        v.np()
    m = v.np_masked()
    assert m.mask.tolist() == [False, True, False] and m.data[2] == 3
    assert RVector("dbl", [1.0, None]).np().tolist()[0] == 1.0
    assert RVector("lgl", [True]).np().dtype == np.bool_
    with pytest.raises(ValueError):
        RVector("int", [1, 2]).scalar()


def test_p_vector_from_dtypes():
    assert p_vector_from(np.array([1, 2], dtype=np.int32))["type"] == "int"
    assert p_vector_from(np.array([1.5]))["type"] == "dbl"
    assert p_vector_from(np.array([True]))["type"] == "lgl"
    assert p_vector_from(np.array(["a", "b"]))["type"] == "chr"
    assert p_vector_from(["a", None])["values"] == ["a", None]
    assert p_vector_from([1, None, 3])["values"] == [1, None, 3]
    assert p_vector_from([True, False])["type"] == "lgl"


def test_p_vector_na_nan_inf_and_masked():
    d = p_vector("dbl", [float("nan"), float("inf"), -float("inf"), None, 1.0])
    assert d["values"] == ["NaN", "Inf", "-Inf", None, 1.0]
    m = np.ma.masked_array([1, 2, 3], mask=[False, True, False])
    assert p_vector_from(m)["values"] == [1, None, 3]
    assert p_vector("int", [float("nan"), 2.0])["values"] == [None, 2]
    with pytest.raises(ValueError):
        p_vector("int", [1.5])


def test_p_vector_names_rules():
    assert p_vector("int", [], names=[])["names"] is None
    assert p_vector("int", [1], names=["a"])["names"] == ["a"]
    with pytest.raises(ValueError):
        p_vector("int", [1, 2], names=["a"])


def test_p_vector_from_cannot_infer_empty():
    with pytest.raises(ValueError):
        p_vector_from([])


def test_p_iranges_granges_seqinfo_shapes():
    ir = p_iranges(np.array([1, 5]), [3, 0], names=["a", "b"], mcols=[("m", p_vector("int", [1, 2]))])
    assert ir["kind"] == "IRanges" and ir["mcols"][0][0] == "m"
    si = p_seqinfo(["a", "b"], [10, None], [False, None], ["hg", "hg"])
    assert si["genome"] == ["hg", "hg"] and si["isCircular"] == [False, None]
    si = p_seqinfo(["a", "b"], [10, None], None, None)
    assert si["isCircular"] == [None, None] and si["seqlengths"] == [10, None]
    gr = p_granges(["a"], [1], [2], ["+"], seqinfo=si)
    assert gr["kind"] == "GRanges" and gr["names"] is None
    with pytest.raises(ValueError):
        p_granges(["a"], [1], [2], ["+"])


def test_p_factor_bases_and_integerlist_shift():
    assert p_factor(["a", "b"], [0, 1, -1], base=0)["codes"] == [1, 2, None]
    assert p_factor(["a", "b"], [1, None, 2])["codes"] == [1, None, 2]
    assert p_integerlist([[0, 1], [], [2]], names=None, shift=1)["elements"] == [[1, 2], [], [3]]


def test_p_hits_and_dataframe_and_null():
    assert p_hits([1, 2], [3, 3], 2, 3) == {"kind": "Hits", "from": [1, 2], "to": [3, 3], "nLnode": 2, "nRnode": 3}
    assert p_dataframe(2, [], None)["rownames"] is None
    assert p_null() == {"kind": "NULL"}


def test_p_rle_runs_are_maximal():
    assert p_rle(np.array([0, 0, 1, 1, 1, 0]))["values"] == [0, 1, 0]
    assert p_rle(np.array([0, 0, 1, 1, 1, 0]))["lengths"] == [2, 3, 1]
    assert p_rle([0, 0, 1], [2, 3, 1])["lengths"] == [5, 1]  # adjacent equal runs merge
    assert p_rle([1, 2], [0, 4])["values"] == [2]  # zero-length runs vanish
    assert p_rle(np.array([1.0, np.nan, np.nan]))["lengths"] == [1, 2]
    na = p_rle(np.ma.masked_array([1, 2, 2], mask=[False, True, True]))
    assert na["values"] == [1, None] and na["lengths"] == [1, 2]
    assert p_rle([], [], type="int") == {"kind": "Rle", "type": "int", "values": [], "lengths": []}


def test_p_rlelist():
    r = p_rlelist(["a", "b"], [np.array([0, 0]), ([1, 2], [1, 1])])
    assert r["names"] == ["a", "b"] and r["elements"][0] == {"type": "int", "values": [0], "lengths": [2]}
    assert p_rlelist([], [])["names"] is None


class Ad(AdapterBase):
    port = "t"
    unsupported_ops = {"merge": "no"}

    def op_findOverlaps(self, *, query, subject=None, maxgap=None, **kw):
        return ("fo", query, subject, maxgap, kw)

    def op_start(self, *, x):
        return ("start", x)

    def op_reduce(self, *, x):
        return ("reduce", x)

    def op_c(self, *dots, x, ignore_mcols=None):
        return ("c", x, dots, ignore_mcols)

    def op_IRanges(self, *, start=None, _mcols=None):
        return ("ir", start, _mcols)


def test_call_binds_and_evaluates_innermost_first():
    ad = Ad()
    env = {"a": "A", "b": "B"}
    call = {"op": "findOverlaps", "args": [{"$ref": "a"}, {"$call": {"op": "start", "args": [{"$ref": "b"}], "kwargs": {}}}, 2],
            "kwargs": {"type": "within", "ignore.strand": True}}
    out = ad.call(call, env)
    assert out[1] == "A" and out[2] == ("start", "B") and out[3].scalar() == 2
    assert out[4]["type"].scalar() == "within" and out[4]["ignore_strand"].scalar() is True


def test_call_dots_mcols_and_errors():
    ad = Ad()
    out = ad.call({"op": "c", "args": [{"$ref": "x"}, {"$ref": "y"}], "kwargs": {"ignore.mcols": True}}, {"x": 1, "y": 2})
    assert out[:3] == ("c", 1, (2,))
    out = ad.call({"op": "IRanges", "args": [], "kwargs": {"start": {"int": [1]}, "score": {"int": [5]}}}, {})
    assert out[2][0][0] == "score"
    with pytest.raises(Unsupported):
        ad.call({"op": "merge", "args": [], "kwargs": {}}, {})
    with pytest.raises(BindError):
        ad.call({"op": "start", "args": [{"$ref": "x"}], "kwargs": {"bogus": 1}}, {"x": 1})
    with pytest.raises(UnknownArgument):  # legal R argument the adapter signature rejects
        ad.call({"op": "reduce", "args": [{"$ref": "x"}], "kwargs": {"with.revmap": True}}, {"x": 1})
    with pytest.raises(AttributeError):  # neither op_ method nor unsupported_ops entry
        ad.call({"op": "end", "args": [{"$ref": "x"}], "kwargs": {}}, {"x": 1})


def test_base_projects_rvector_and_none():
    ad = Ad()
    assert ad._project(None) == {"kind": "NULL"}
    assert ad._project(RVector("int", [1, None], ["a", "b"])) == {"kind": "vector", "type": "int", "values": [1, None], "names": ["a", "b"]}
