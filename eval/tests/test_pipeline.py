"""End to end on a slice of the pool: oracle -> r_echo / mutant through run_port.py -> comparator -> report."""

import json
import shutil
import subprocess
import sys
from collections import Counter

import pytest

from bioc_eval import oracle, report
from bioc_eval.common import EVAL, iter_jsonl, load_pool, write_jsonl
from bioc_eval.compare import score

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="needs Rscript")


@pytest.fixture(scope="module")
def slice_build(tmp_path_factory):
    out = tmp_path_factory.mktemp("slice")
    pool = load_pool()
    ids = {c["id"] for i, c in enumerate(pool) if i % 11 == 0}
    meta = oracle.build(workers=2, out_dir=out, subset=ids)
    assert meta["cases"] == len(ids)
    return out


def run_adapter(name, out):
    raw = out / f"{name}.jsonl"
    subprocess.run([sys.executable, str(EVAL / "runner" / "run_port.py"), "--port", name, "--resolved", str(out / "resolved.jsonl"),
                    "--expected", str(out / "expected.jsonl"), "--out", str(raw)], check=True, capture_output=True)
    return raw


def scored(out, raw):
    exp = {r["id"]: r for r in iter_jsonl(out / "expected.jsonl") if r["status"] == "ok"}
    got = {r["id"]: r for r in iter_jsonl(raw)}
    return [score(exp[i], got.get(i)) for i in sorted(exp)]


def test_r_echo_scores_100_percent(slice_build):
    rows = scored(slice_build, run_adapter("r_echo", slice_build))
    assert len(rows) > 500 and {r["outcome"] for r in rows} == {"pass"}


def test_mutant_outcomes_equal_predictions(slice_build):
    raw = run_adapter("mutant", slice_build)
    pred = {r["id"]: r for r in iter_jsonl(slice_build / "mutant.predicted.jsonl")}
    rows = scored(slice_build, raw)
    assert [r["id"] for r in rows if pred[r["id"]]["predicted"] != r["outcome"]] == []
    assert len({p["perturbation"] for p in pred.values() if p["perturbation"]}) >= 7
    assert {"mismatch", "adapter_error", "unsupported", "port_error_unexpected"} <= {r["outcome"] for r in rows}


def test_resolved_has_no_expected_values_and_only_projections(slice_build):
    for rec in iter_jsonl(slice_build / "resolved.jsonl"):
        assert set(rec) == {"id", "source", "op", "hash", "inputs", "call"}
        assert all(p["kind"] in {"vector", "IRanges", "GRanges", "Seqinfo", "DataFrame", "Factor", "IntegerList"}
                   for p in rec["inputs"].values())


def test_runner_resume_skips_finished_ids(slice_build):
    raw = slice_build / "r_echo.jsonl"
    n = len(raw.read_text().splitlines())
    r = subprocess.run([sys.executable, str(EVAL / "runner" / "run_port.py"), "--port", "r_echo", "--resolved", str(slice_build / "resolved.jsonl"),
                        "--expected", str(slice_build / "expected.jsonl"), "--out", str(raw), "--resume"], capture_output=True, text=True)
    assert r.returncode == 0 and len(raw.read_text().splitlines()) == n and f"ran 0 cases ({n} resumed)" in r.stderr


def test_report_reconciles_and_renders_port_sections(slice_build, tmp_path):
    rows = scored(slice_build, run_adapter("mutant", slice_build))
    (tmp_path / "scored").mkdir()
    write_jsonl(tmp_path / "scored" / "mutant.jsonl", rows)
    write_jsonl(tmp_path / "scored" / "r_echo.jsonl", scored(slice_build, run_adapter("r_echo", slice_build)))
    expected = list(iter_jsonl(slice_build / "expected.jsonl"))
    summary = {}
    for suite in ("e1", "e2", "e4"):
        n = sum(1 for e in expected if e["source"].startswith(suite))
        summary[suite] = {"read": n, "converted": n, "unconvertible": 0}
    (tmp_path / "conversion-summary.json").write_text(json.dumps(summary))
    (tmp_path / "unconvertible.jsonl").write_text("")
    registry = {"mutant": {"tag": None}, "r_echo": {"selftest": True}}  # mutant plays a real port here
    text, data = report.build(results=tmp_path, build_dir=slice_build, cases=tmp_path, registry=registry)
    assert "## Headline: pass rate per port x source" in text and "## mutant: family x source" in text
    assert "## Self-tests" in text and "| r_echo |" in text
    assert set(data["matrix"]["mutant"]) >= {"overlaps", "intra", "construct"}
    total = sum(sum(c.values()) for srcs in data["matrix"]["mutant"].values() for c in srcs.values())
    assert total == len(rows)
    # reconciliation is enforced: a wrong 'read' count must fail loudly
    summary["e1"]["read"] += 1
    (tmp_path / "conversion-summary.json").write_text(json.dumps(summary))
    with pytest.raises(AssertionError):
        report.build(results=tmp_path, build_dir=slice_build, cases=tmp_path, registry=registry)
