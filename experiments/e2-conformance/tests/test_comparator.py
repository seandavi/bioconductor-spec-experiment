"""The comparator must reject differences, or the whole suite proves nothing."""

from __future__ import annotations

from ranges_e2.conformance.codec import compare, enc_vec, tree_equal


def test_rejects_value_and_type_differences() -> None:
    assert not tree_equal(enc_vec("int", [1, 2]), enc_vec("int", [1, 3]))
    assert not tree_equal(enc_vec("int", [1]), enc_vec("dbl", [1.0]))
    assert not tree_equal(enc_vec("int", [1, None]), enc_vec("int", [1, 0]))
    assert not tree_equal(enc_vec("lgl", [True]), enc_vec("int", [1]))
    assert not tree_equal({"a": 1, "b": 2}, {"b": 2, "a": 1})  # key order is significant


def test_double_tolerance() -> None:
    assert tree_equal(enc_vec("dbl", [1.0]), enc_vec("dbl", [1.0 + 1e-12]))
    assert not tree_equal(enc_vec("dbl", [1.0]), enc_vec("dbl", [1.001]))


def test_hits_query_grouped_mode() -> None:
    def hits(frm: list[int], to: list[int]) -> dict[str, object]:
        return {"$type": "Hits", "from": enc_vec("int", frm), "to": enc_vec("int", to), "nLnode": 3, "nRnode": 3}

    assert compare(hits([1, 1, 2], [3, 1, 2]), hits([1, 1, 2], [1, 3, 2]), "hits-query-grouped")
    assert not compare(hits([1, 1, 2], [3, 1, 2]), hits([1, 1, 2], [1, 3, 2]), "exact")
    assert not compare(hits([1, 2, 1], [1, 2, 3]), hits([1, 1, 2], [1, 3, 2]), "hits-query-grouped")
