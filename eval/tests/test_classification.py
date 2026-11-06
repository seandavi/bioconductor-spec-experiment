"""Runner classification: traceback-based split between port errors and adapter errors."""

import time

import fakeport
from adapter_base import AdapterBase, RVector, Unsupported, UnknownArgument
from binding import BindError
from run_port import classify, enters_port, run_one


class Probe(AdapterBase):
    port = "probe"
    port_packages = ("fakeport",)
    unsupported_ops = {"merge": "nope"}

    def __init__(self, behavior):
        self.behavior = behavior

    def run_case(self, case, stage):
        stage("build", lambda: None)
        stage("call", self.behavior)
        return {"result": {"kind": "NULL"}}


CASE = {"id": "c", "inputs": {}, "call": {"op": "start", "args": [], "kwargs": {}}}


def run(behavior, timeout=5.0):
    return run_one(Probe(behavior), CASE, timeout)


def test_value():
    r = run(lambda: None)
    assert r["status"] == "value" and r["result"] == {"kind": "NULL"}


def test_exception_through_port_frames_is_port_error():
    r = run(lambda: fakeport.boom("x"))
    assert (r["status"], r["phase"], r["error"]["origin"]) == ("port_error", "call", "port")


def test_exception_in_adapter_code_is_adapter_error():
    def adapter_bug():
        return len(5)  # TypeError in translation

    r = run(adapter_bug)
    assert r["status"] == "adapter_error" and r["error"]["type"] == "TypeError"


def test_unsupported_and_unknown_argument():
    def unsupported():
        raise Unsupported("subscript:negative")

    def unknown():
        raise UnknownArgument("nthread")

    assert run(unsupported)["status"] == "unsupported"
    r = run(unknown)
    assert r["status"] == "unsupported" and r["reason"] == "argument:nthread"


def test_bind_error_is_a_port_error_with_origin_bind():
    def bad():
        raise BindError("unused argument (x)")

    r = run(bad)
    assert (r["status"], r["error"]["origin"]) == ("port_error", "bind")


def test_timeout_is_a_port_error_even_if_port_code_swallows_exceptions():
    def hang():
        try:
            time.sleep(5)
        except Exception:  # a port's broad handler must not swallow the timeout
            pass

    r = run(hang, timeout=0.2)
    assert r["status"] == "port_error" and r["error"]["origin"] == "timeout"


def test_build_phase_is_recorded():
    class B(Probe):
        def run_case(self, case, stage):
            stage("build", fakeport.boom, "bad input")
            return {"result": {}}

    r = run_one(B(None), CASE, 5.0)
    assert (r["status"], r["phase"]) == ("port_error", "build")


def test_enters_port_looks_at_module_names():
    try:
        fakeport.boom()
    except RuntimeError as e:
        assert enters_port(e.__traceback__, ("fakeport",))
        assert not enters_port(e.__traceback__, ("other",))
        assert classify(e, "call", ("fakeport",))["status"] == "port_error"


def test_projection_failures_are_phase_project_and_never_pass_against_an_r_error():
    """Op succeeds, then a native getter fails (or hangs) during projection; R errored on the same case."""
    from bioc_eval.compare import score

    class P(Probe):
        def __init__(self, project_fn):
            self.project_fn = project_fn

        def run_case(self, case, stage):
            stage("build", lambda: None)
            value = stage("call", lambda: object())
            return {"result": stage("project", self.project_fn, value)}

    exp = {"id": "c", "source": "heldout", "op": "IRanges", "status": "ok",
           "result": {"kind": "error", "message": "each range must have a non-negative width"}}

    boom = run_one(P(lambda v: fakeport.boom("getter failed")), CASE, 5.0)
    assert (boom["status"], boom["phase"]) == ("port_error", "project")
    assert score(exp, boom)["outcome"] == "port_error_unexpected"

    def hang(v):
        time.sleep(5)

    slow = run_one(P(hang), CASE, 0.2)
    assert (slow["phase"], slow["error"]["origin"]) == ("project", "timeout")
    assert score(exp, slow)["outcome"] == "port_error_unexpected"

    bug = run_one(P(lambda v: len(5)), CASE, 5.0)  # adapter-side projection bug
    assert score(exp, bug)["outcome"] == "adapter_error"
