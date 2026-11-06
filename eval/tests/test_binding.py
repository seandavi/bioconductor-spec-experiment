import pytest

from binding import BindError, bind
from bioc_eval.common import load_ops

OPS = load_ops()


def test_positional_to_names_and_named():
    b, pos, named = bind(OPS["findOverlaps"], ["Q", "S", 3], {"type": "within", "ignore.strand": True})
    assert b == {"query": "Q", "subject": "S", "maxgap": 3, "type": "within", "ignore.strand": True}
    assert pos == [] and named == []


def test_objects_dots_after_positional():
    b, pos, _ = bind(OPS["c"], ["X", "Y", "Z"], {"ignore.mcols": True})
    assert b == {"x": "X", "ignore.mcols": True} and pos == ["Y", "Z"]


def test_order_has_no_positional_formals():
    b, pos, _ = bind(OPS["order"], ["X", "Y"], {"decreasing": True})
    assert b == {"decreasing": True} and pos == ["X", "Y"]


def test_mcols_dots_collect_unknown_named():
    b, pos, named = bind(OPS["IRanges"], [1, 2], {"width": 3, "score": "S"})
    assert b == {"start": 1, "end": 2, "width": 3}
    assert named == [("score", "S")]


def test_unused_argument_is_a_bind_error():
    with pytest.raises(BindError):
        bind(OPS["start"], ["X"], {"bogus": 1})
    with pytest.raises(BindError):
        bind(OPS["start"], ["X", "extra"], {})


def test_duplicate_match_is_a_bind_error():
    with pytest.raises(BindError):
        bind(OPS["shift"], ["X", 1], {"shift": 2})
