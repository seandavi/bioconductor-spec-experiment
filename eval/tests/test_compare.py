"""Every comparator row and rule (spec/COMPARATOR.md)."""

from bioc_eval.compare import diffs, score


def vec(vals, type_="int", names=None):
    return {"kind": "vector", "type": type_, "values": vals, "names": names}


def hits(f, t, nl=3, nr=3):
    return {"kind": "Hits", "from": f, "to": t, "nLnode": nl, "nRnode": nr}


def exp(result, **extra):
    return {"id": "x", "source": "heldout", "op": "start", "status": "ok", "result": result, **extra}


def value(result):
    return {"status": "value", "result": result}


ERR = {"kind": "error", "message": "boom"}


def outcome(e, rec):
    return score(e, rec)["outcome"]


# ---- table rows
def test_row1_unsupported_even_when_r_errors():
    assert outcome(exp(ERR), {"status": "unsupported", "reason": "x"}) == "unsupported"


def test_row2_adapter_error_never_passes():
    assert outcome(exp(ERR), {"status": "adapter_error", "phase": "call", "error": {}}) == "adapter_error"


def test_row2b_build_port_error_is_unexpected_even_if_r_errors():
    assert outcome(exp(ERR), {"status": "port_error", "phase": "build", "error": {}}) == "port_error_unexpected"


def test_row3_both_error_passes_R1():
    r = score(exp(ERR), {"status": "port_error", "phase": "call", "error": {"message": "different"}})
    assert (r["outcome"], r["rule"]) == ("pass", "R1")


def test_row2c_projection_failure_never_earns_R1():
    # the operation returned, then reading the result back failed inside the port: R errored, so this is no "both error"
    r = score(exp(ERR), {"status": "port_error", "phase": "project", "error": {"origin": "port"}})
    assert (r["outcome"], r["rule"]) == ("port_error_unexpected", "row2c")
    r = score(exp(ERR), {"status": "port_error", "phase": "project", "error": {"origin": "timeout"}})
    assert r["outcome"] == "port_error_unexpected"


def test_row4b_timeout_in_call_does_not_match_an_r_error():
    r = score(exp(ERR), {"status": "port_error", "phase": "call", "error": {"origin": "timeout"}})
    assert r["outcome"] == "port_error_unexpected"
    r = score(exp(ERR), {"status": "port_error", "phase": "call", "error": {"origin": "bind"}})
    assert (r["outcome"], r["rule"]) == ("pass", "R1")


def test_row4_r_error_port_value_mismatch():
    r = score(exp(ERR), value(vec([1])))
    assert r["outcome"] == "mismatch" and r["detail"]["kinds"] == ["r_error_port_value"]


def test_row5_port_error_where_r_succeeded():
    assert outcome(exp(vec([1])), {"status": "port_error", "phase": "call", "error": {"origin": "timeout"}}) == "port_error_unexpected"


def test_missing_record_is_unsupported():
    assert outcome(exp(vec([1])), None) == "unsupported"


def test_invalid_projection_is_adapter_error():
    assert outcome(exp(vec([1])), value({"kind": "vector", "type": "int", "values": [1.5], "names": None})) == "adapter_error"


# ---- exact
def test_exact_pass_and_dbl_equals():
    assert outcome(exp(vec([1.0, 0.0], "dbl")), value(vec([1.0, -0.0], "dbl"))) == "pass"
    assert outcome(exp(vec(["NaN", "Inf", None], "dbl")), value(vec(["NaN", "Inf", None], "dbl"))) == "pass"


def test_null_equals_only_null():
    assert outcome(exp(vec([None], "lgl")), value(vec([False], "lgl"))) == "mismatch"


def test_bool_is_not_int():
    assert outcome(exp(vec([1], "int")), value(vec([True], "lgl"))) == "mismatch"


def test_names_and_type_are_not_relaxed():
    r = score(exp(vec([1], "int", ["a"])), value(vec([1], "int", None)))
    assert r["outcome"] == "mismatch" and r["detail"]["kinds"] == ["names"]
    r = score(exp(vec([1], "int")), value(vec([1.0], "dbl")))
    assert r["outcome"] == "mismatch" and "type" in r["detail"]["kinds"]


def test_length_and_value_and_order_kinds():
    assert score(exp(vec([1, 2])), value(vec([1])))["detail"]["kinds"] == ["length"]
    assert score(exp(vec([1, 2])), value(vec([1, 3])))["detail"]["kinds"] == ["value"]
    assert score(exp(vec([1, 2])), value(vec([2, 1])))["detail"]["kinds"] == ["order"]


def test_float_ulp_only_when_every_difference_is():
    a = 0.1 + 0.2
    r = score(exp(vec([a], "dbl")), value(vec([0.3], "dbl")))
    assert r["outcome"] == "mismatch" and r["detail"]["kinds"] == ["float_ulp"]
    r = score(exp(vec([a, 5.0], "dbl")), value(vec([0.3, 6.0], "dbl")))
    assert r["detail"]["kinds"] == ["value"]


def test_kind_difference():
    r = score(exp(vec([1])), value({"kind": "NULL"}))
    assert r["detail"]["kinds"] == ["kind"]


def test_diffs_reports_first_path():
    ds = diffs({"a": [1, {"b": 2}]}, {"a": [1, {"b": 3}]})
    assert ds[0]["path"] == "$.a[1].b"


# ---- R2: Hits as a set within a query, query order exact
def test_R2_subject_order_within_query_free():
    e = exp(hits([1, 1, 2], [1, 3, 2]), sorted_by_query=True)
    assert outcome(e, value(hits([1, 1, 2], [3, 1, 2]))) == "pass"


def test_R2_query_order_still_exact():
    e = exp(hits([1, 2], [1, 2]), sorted_by_query=True)
    assert outcome(e, value(hits([2, 1], [2, 1]))) == "mismatch"


def test_R2_missing_and_extra_pairs_and_dims():
    e = exp(hits([1, 1], [1, 2]), sorted_by_query=True)
    assert outcome(e, value(hits([1], [1]))) == "mismatch"
    assert outcome(e, value(hits([1, 1], [1, 3]))) == "mismatch"
    assert outcome(e, value(hits([1, 1], [1, 2], nr=4))) == "mismatch"


def test_R2_not_applied_when_r_class_is_not_sorted_by_query():
    e = exp(hits([1, 1], [1, 2]), sorted_by_query=False)
    assert outcome(e, value(hits([1, 1], [2, 1]))) == "mismatch"


# ---- R3: select="arbitrary"
U = hits([1, 1, 2], [1, 3, 2], 3, 3)


def test_R3_any_universe_hit_passes():
    e = exp(vec([1, 2, None]), universe=U)
    assert outcome(e, value(vec([3, 2, None]))) == "pass"


def test_R3_outside_universe_and_wrong_null_fail():
    e = exp(vec([1, 2, None]), universe=U)
    assert score(e, value(vec([2, 2, None])))["detail"]["kinds"] == ["not_in_universe"]
    assert outcome(e, value(vec([1, 2, 1]))) == "mismatch"
    assert outcome(e, value(vec([None, 2, None]))) == "mismatch"


def test_R3_rs_own_pick_is_allowed():
    e = exp(vec([1]), universe=hits([1], [2], 1, 3))
    assert outcome(e, value(vec([1]))) == "pass"
    assert outcome(e, value(vec([3]))) == "mismatch"


def test_R3_requires_names_type_and_length():
    e = exp(vec([1], "int", ["q"]), universe=hits([1], [1], 1, 1))
    assert outcome(e, value(vec([1]))) == "mismatch"
    assert outcome(exp(vec([1]), universe=hits([1], [1], 1, 1)), value(vec([1, 1]))) == "mismatch"
