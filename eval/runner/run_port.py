"""Port runner: executes resolved cases through one adapter inside the port's own venv.

    uv run --frozen --project <port dir> python eval/runner/run_port.py --port NAME \
        --resolved build/resolved.jsonl --out build/results/NAME.jsonl [--resume] [--timeout 30]

stdlib + numpy + the port only. Never reads expected values (self-test adapters get them via --expected).
Output: one JSON record per case, flushed:
  {id, status: value|port_error|unsupported|adapter_error, phase, result?, error?, reason?, warnings, ms}
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import json
import signal
import sys
import time
import traceback
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVAL = HERE.parent
sys.path.insert(0, str(HERE))

from adapter_base import AdapterBase, Unsupported, UnknownArgument  # noqa: E402
from binding import BindError  # noqa: E402


class CaseTimeout(BaseException):
    """Raised by SIGALRM. BaseException so port code with `except Exception` cannot swallow it."""


def load_adapter(name: str, ctx: dict) -> AdapterBase:
    sys.path.insert(0, str(EVAL / "adapters"))
    if ctx.get("selftest_path"):
        sys.path.insert(0, str(EVAL / "selftest"))
    path = EVAL / "adapters" / f"{name}.py"
    if not path.exists():
        raise SystemExit(f"no adapter for port {name!r}: expected {path}")
    spec = importlib.util.spec_from_file_location(f"adapter_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    adapter = mod.Adapter()
    if hasattr(adapter, "setup"):
        adapter.setup(ctx)
    return adapter


def enters_port(tb, packages: tuple[str, ...]) -> bool:
    """Does the traceback pass through a frame whose module is in one of `packages`?"""
    while tb is not None:
        mod = tb.tb_frame.f_globals.get("__name__", "")
        if any(mod == p or mod.startswith(p + ".") for p in packages):
            return True
        tb = tb.tb_next
    return False


def classify(exc: BaseException, phase: str, packages: tuple[str, ...]) -> dict:
    """Mechanical split (ADAPTERS.md, Classification)."""
    if isinstance(exc, Unsupported):
        return {"status": "unsupported", "reason": str(exc)[:300]}
    if isinstance(exc, UnknownArgument):
        return {"status": "unsupported", "reason": f"argument:{exc.name}"}
    err = {"type": type(exc).__name__, "message": str(exc)[:500]}
    if isinstance(exc, CaseTimeout):
        return {"status": "port_error", "phase": phase, "error": {**err, "origin": "timeout"}}
    if isinstance(exc, BindError):
        return {"status": "port_error", "phase": phase, "error": {**err, "origin": "bind"}}
    if enters_port(exc.__traceback__, packages):
        return {"status": "port_error", "phase": phase, "error": {**err, "origin": "port"}}
    tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__)[-4:])
    return {"status": "adapter_error", "phase": phase, "error": {**err, "traceback": tb[-800:]}}


def run_one(adapter: AdapterBase, case: dict, timeout: float) -> dict:
    phase = ["build"]

    def stage(name, fn, *a):
        phase[0] = name
        return fn(*a)

    def on_alarm(signum, frame):
        raise CaseTimeout(f"case exceeded {timeout} s")

    t0 = time.perf_counter()
    rec: dict
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        signal.signal(signal.SIGALRM, on_alarm)
        signal.setitimer(signal.ITIMER_REAL, timeout)
        try:
            out = adapter.run_case(case, stage)
            rec = {"status": "value", **out}
        except (Exception, CaseTimeout) as e:  # noqa: BLE001 - classification is the point
            rec = classify(e, phase[0], adapter.port_packages)
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
    if "warnings" not in rec:
        rec["warnings"] = sorted({f"{w.category.__name__}: {w.message}"[:200] for w in caught})[:10]
    rec["id"] = case["id"]
    rec["ms"] = round((time.perf_counter() - t0) * 1000, 2)
    return rec


def read_done(path: Path) -> set[str]:
    if not path.exists():
        return set()
    done = set()
    for line in path.read_text().splitlines():
        try:
            done.add(json.loads(line)["id"])
        except (json.JSONDecodeError, KeyError):
            pass  # truncated last line after a kill
    return done


def roundtrip(adapter: AdapterBase, resolved: list[dict]) -> int:
    """project(build(p)) == p for every resolved input; prints failures per kind."""
    seen: dict[str, dict] = {}
    for case in resolved:
        for name, p in case["inputs"].items():
            seen.setdefault(json.dumps(p, sort_keys=True), p)
    bad: dict[str, list] = {}
    n = 0
    for key, p in seen.items():
        n += 1
        try:
            back = adapter._project(adapter.build_input(p))
            if json.dumps(back, sort_keys=True) != key:
                bad.setdefault(p["kind"], []).append(("differs", key[:160], json.dumps(back, sort_keys=True)[:160]))
        except Unsupported as e:
            bad.setdefault(p["kind"], []).append(("unsupported", str(e), ""))
        except Exception as e:  # noqa: BLE001
            where = "port" if enters_port(e.__traceback__, adapter.port_packages) else "adapter"
            bad.setdefault(p["kind"], []).append((f"{where}:{type(e).__name__}", str(e)[:120], ""))
    print(f"roundtrip: {n} distinct inputs")
    for kind, items in sorted(bad.items()):
        print(f"  {kind}: {len(items)} failing")
        for item in items[:3]:
            print("     ", item)
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", required=True)
    ap.add_argument("--resolved", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--expected", type=Path, help="self-test adapters only")
    ap.add_argument("--timeout", type=float, default=30.0)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--roundtrip", action="store_true")
    ns = ap.parse_args()

    ctx = {"expected_path": str(ns.expected) if ns.expected else None, "selftest_path": bool(ns.expected),
           "out_path": str(ns.out) if ns.out else None}
    adapter = load_adapter(ns.port, ctx)
    cases = [json.loads(l) for l in ns.resolved.read_text().splitlines() if l.strip()]
    if ns.roundtrip:
        return roundtrip(adapter, cases)
    if ns.limit:
        cases = cases[: ns.limit]
    ns.out.parent.mkdir(parents=True, exist_ok=True)
    done = read_done(ns.out) if ns.resume else set()
    mode = "a" if ns.resume else "w"
    n = 0
    with ns.out.open(mode) as fh:
        for case in cases:
            if case["id"] in done:
                continue
            fh.write(json.dumps(run_one(adapter, case, ns.timeout), ensure_ascii=False, separators=(",", ":")) + "\n")
            fh.flush()
            n += 1
    if hasattr(adapter, "finish"):
        adapter.finish()
    print(f"{ns.port}: ran {n} cases ({len(done)} resumed) -> {ns.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
