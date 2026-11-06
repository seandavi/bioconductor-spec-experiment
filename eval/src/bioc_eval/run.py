"""`just run PORT`: verify the frozen port, run it in its own venv, score against expected.jsonl."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

from .common import BUILD, EVAL, RESULTS, WORKTREES, iter_jsonl, ports, write_jsonl
from .compare import OUTCOMES, score


def git(repo: Path, *args: str) -> str:
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, check=True).stdout.strip()


def port_state(spec: dict) -> dict:
    """HEAD sha and cleanliness of the port's worktree (self-tests and biocpy have no frozen sha)."""
    if "sha" not in spec:
        return {"head": None, "clean": True}
    repo = WORKTREES / spec["project"].split("/")[0]
    return {"head": git(repo, "rev-parse", "--short", "HEAD"), "clean": git(repo, "status", "--porcelain") == ""}


def check_frozen(name: str, spec: dict, when: str) -> dict:
    state = port_state(spec)
    if "sha" in spec and (not state["head"].startswith(spec["sha"]) or not state["clean"]):
        raise SystemExit(f"{when}: port {name} is not frozen: HEAD {state['head']} (want {spec['sha']}), "
                         f"clean={state['clean']}")
    return state


def project_dir(spec: dict) -> Path:
    """`eval` = this tools project; `eval/...` is relative to the repo root; else relative to WORKTREES."""
    proj = spec["project"]
    if proj == "eval":
        return EVAL
    if proj.startswith("eval/"):
        return EVAL.parent / proj
    return WORKTREES / proj


def run_command(name: str, spec: dict, out: Path, extra: list[str]) -> list[str]:
    cmd = ["uv", "run", "--frozen", "--project", str(project_dir(spec)), "python", str(EVAL / "runner" / "run_port.py"),
           "--port", name, "--resolved", str(BUILD / "resolved.jsonl"), "--out", str(out)]
    if spec.get("selftest"):
        cmd += ["--expected", str(BUILD / "expected.jsonl")]
    return cmd + extra


def adapter_digest(name: str) -> str | None:
    p = EVAL / "adapters" / f"{name}.py"
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def lock_digest(spec: dict) -> str | None:
    lock = project_dir(spec) / "uv.lock"
    return hashlib.sha256(lock.read_bytes()).hexdigest() if lock.exists() else None


def score_file(raw: Path, scored: Path) -> Counter:
    expected = {r["id"]: r for r in iter_jsonl(BUILD / "expected.jsonl") if r["status"] == "ok"}
    got = {r["id"]: r for r in iter_jsonl(raw)}
    rows = [score(e, got.get(i)) for i, e in sorted(expected.items())]
    write_jsonl(scored, rows)
    return Counter(r["outcome"] for r in rows)


def adapter_check(raw: Path) -> None:
    """Only adapter_error / unsupported (and what the port raised): no R values are read."""
    recs = list(iter_jsonl(raw))
    by_status = Counter(r["status"] for r in recs)
    print(f"{len(recs)} cases: {dict(by_status)}")
    for status, key in (("adapter_error", lambda r: (r["phase"], r["error"]["type"], r["error"]["message"][:100])),
                        ("unsupported", lambda r: (r.get("reason", "")[:120],)),
                        ("port_error", lambda r: (r["phase"], r["error"]["type"], r["error"]["message"][:80]))):
        c = Counter(key(r) for r in recs if r["status"] == status)
        if c:
            print(f"\n{status} ({sum(c.values())}):")
            for k, n in c.most_common(25):
                print(f"  {n:5d}  {' | '.join(map(str, k))}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("port")
    ap.add_argument("--roundtrip", action="store_true")
    ap.add_argument("--adapter-check", action="store_true")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--timeout", type=float, default=30.0)
    ap.add_argument("--resume", action="store_true")
    ns = ap.parse_args()
    reg = ports()
    if ns.port not in reg:
        raise SystemExit(f"unknown port {ns.port!r}; known: {', '.join(reg)}")
    spec = reg[ns.port]
    raw = BUILD / "results" / f"{ns.port}.jsonl"
    extra = ["--timeout", str(ns.timeout)]
    if ns.limit:
        extra += ["--limit", str(ns.limit)]
    if ns.resume:
        extra += ["--resume"]
    if ns.roundtrip:
        cmd = run_command(ns.port, spec, raw, extra + ["--roundtrip"])
        return subprocess.run(cmd, cwd=EVAL).returncode
    before = check_frozen(ns.port, spec, "before")
    t0 = time.time()
    rc = subprocess.run(run_command(ns.port, spec, raw, extra), cwd=EVAL).returncode
    seconds = round(time.time() - t0, 1)
    after = check_frozen(ns.port, spec, "after")
    if rc != 0:
        print(f"runner exited {rc}", file=sys.stderr)
        return rc
    if ns.adapter_check:
        adapter_check(raw)
        return 0
    scored = RESULTS / "scored" / f"{ns.port}.jsonl"
    counts = score_file(raw, scored)
    oracle_meta = json.loads((BUILD / "oracle-meta.json").read_text())
    meta = {"port": ns.port, "tag": spec.get("tag"), "sha_before": before["head"], "sha_after": after["head"],
            "clean_before": before["clean"], "clean_after": after["clean"], "adapter_sha256": adapter_digest(ns.port), "lock_sha256": lock_digest(spec),
            "run_seconds": seconds, "oracle_versions": oracle_meta["versions"], "outcomes": dict(counts)}
    (RESULTS / "scored" / f"{ns.port}.meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    total = sum(counts.values())
    print(f"{ns.port}: {total} scored in {seconds}s  " + "  ".join(f"{o}={counts.get(o, 0)}" for o in OUTCOMES))
    return 0


if __name__ == "__main__":
    sys.exit(main())
