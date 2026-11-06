"""`just oracle`: run R (normative) over the whole pool -> build/resolved.jsonl, build/expected.jsonl.

R runs in chunks (parallel Rscript processes). Per-case results are cached by
sha256(R + package versions, oracle.R digest, case hash) in build/oracle-cache.jsonl.
A crashed chunk is resumed after marking the case that killed it `excluded: r_crash`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .common import BUILD, EVAL, RESULTS, case_hash, canonical, iter_jsonl, load_ops, load_pool, write_jsonl

sys.path.insert(0, str(EVAL / "runner"))
from binding import BindError, bind  # noqa: E402

ORACLE_R = EVAL / "oracle" / "oracle.R"
CACHE = BUILD / "oracle-cache.jsonl"
CHUNK = 400


def r_versions() -> str:
    code = ('suppressWarnings(suppressMessages({library(IRanges);library(GenomicRanges);library(jsonlite)}));'
            'v<-function(p) as.character(packageVersion(p));'
            'cat(paste0("R ",getRversion(),"; IRanges ",v("IRanges"),"; GenomicRanges ",v("GenomicRanges"),'
            '"; S4Vectors ",v("S4Vectors"),"; Seqinfo ",v("Seqinfo"),"; BiocGenerics ",v("BiocGenerics"),'
            '"; jsonlite ",v("jsonlite")))')
    return subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=True).stdout.strip()


def oracle_digest() -> str:
    return hashlib.sha256(ORACLE_R.read_bytes()).hexdigest()


def fingerprint(versions: str) -> str:
    return hashlib.sha256(f"{versions}|{oracle_digest()}".encode()).hexdigest()


# ------------------------------------------------------------------ universe (COMPARATOR.md, R3)

def _chr_scalar(v) -> str | None:
    if isinstance(v, str):
        return v
    if isinstance(v, dict) and "chr" in v and len(v["chr"]) == 1:
        return v["chr"][0]
    return None


def universe_call(case: dict, ops: dict) -> dict | None:
    """The same call with select = "all", where rule R3 can apply; else None."""
    call = case["call"]
    if call["op"] not in ("findOverlaps", "nearest"):
        return None
    try:
        bound, extra_pos, extra_named = bind(ops[call["op"]], call["args"], call["kwargs"])
    except BindError:
        return None
    sel = bound.get("select")
    if call["op"] == "findOverlaps":
        applies = "select" in bound and _chr_scalar(sel) == "arbitrary"
    else:
        applies = "select" not in bound or _chr_scalar(sel) == "arbitrary"
    if not applies:
        return None
    kwargs = dict(bound)
    kwargs["select"] = "all"
    return {"op": call["op"], "args": [], "kwargs": kwargs}


# ------------------------------------------------------------------ R execution

def _run_chunk(lines: list[str], workdir: Path, idx: int) -> list[dict]:
    """Run one Rscript over `lines`; resume past crashes. Returns raw R records (in order)."""
    out: list[dict] = []
    remaining = lines
    attempt = 0
    while remaining:
        src = workdir / f"in-{idx}-{attempt}.jsonl"
        dst = workdir / f"out-{idx}-{attempt}.jsonl"
        src.write_text("".join(remaining))
        log = workdir / f"r-{idx}-{attempt}.log"
        with log.open("w") as fh:
            subprocess.run(["Rscript", str(ORACLE_R), str(src), str(dst)], stderr=fh, stdout=fh)
        got = [json.loads(l) for l in dst.read_text().splitlines() if l.strip()] if dst.exists() else []
        out.extend(got)
        done = len(got)
        if done >= len(remaining):
            break
        victim = json.loads(remaining[done])["id"]
        out.append({"id": victim, "status": "excluded", "reason": "r_crash", "detail": "Rscript died", "warnings": []})
        remaining = remaining[done + 1:]
        attempt += 1
    return out


def run_r(cases: list[dict], ops: dict, workers: int) -> dict[str, dict]:
    if not cases:
        return {}
    with tempfile.TemporaryDirectory(prefix="oracle-") as td:
        workdir = Path(td)
        lines = []
        for c in cases:
            rec = {"id": c["id"], "inputs": c["inputs"], "call": c["call"], "universe_call": universe_call(c, ops)}
            lines.append(json.dumps(rec, ensure_ascii=False, separators=(",", ":")) + "\n")
        chunks = [lines[i:i + CHUNK] for i in range(0, len(lines), CHUNK)]
        with ThreadPoolExecutor(max_workers=workers) as ex:
            parts = list(ex.map(lambda ic: _run_chunk(ic[1], workdir, ic[0]), enumerate(chunks)))
        result = {r["id"]: r for part in parts for r in part}
        # keep R stderr for inspection
        logdir = BUILD / "oracle-logs"
        logdir.mkdir(parents=True, exist_ok=True)
        for f in logdir.glob("*.log"):
            f.unlink()
        for f in workdir.glob("r-*.log"):
            (logdir / f.name).write_bytes(f.read_bytes())
    return result


# ------------------------------------------------------------------ features (PROJECTION.md, "Oracle extras")

def _walk_vectors(node):
    if isinstance(node, dict):
        for k in ("int", "dbl", "chr", "lgl"):
            if k in node and isinstance(node[k], list):
                yield node[k]
        for v in node.values():
            yield from _walk_vectors(v)
    elif isinstance(node, list):
        for v in node:
            yield from _walk_vectors(v)


def features_of(case: dict, inputs_proj: dict, errored: bool) -> list[str]:
    f: set[str] = set()
    for p in inputs_proj.values():
        k = p.get("kind")
        if k in ("IRanges", "GRanges"):
            if len(p["start"]) == 0:
                f.add("empty")
            if 0 in p["width"]:
                f.add("zero_width")
            if p["names"] is not None:
                f.add("names")
            if p["mcols"]:
                f.add("mcols")
        if k == "GRanges":
            if len(set(p["strand"])) > 1 or any(s != "*" for s in p["strand"]):
                f.add("strand")
            si = p["seqinfo"]
            if len(si["seqnames"]) > 1:
                f.add("multi_seqlevel")
            if any(v is not None for v in si["seqlengths"]):
                f.add("seqlengths")
            if any(v for v in si["isCircular"]):
                f.add("circular")
        if k == "Seqinfo":
            if len(p["seqnames"]) > 1:
                f.add("multi_seqlevel")
            if any(v is not None for v in p["seqlengths"]):
                f.add("seqlengths")
            if any(v for v in p["isCircular"]):
                f.add("circular")
    for vals in _walk_vectors([case["call"]]):
        if any(v is None or v == "NaN" for v in vals):
            f.add("na_arg")
    if errored:
        f.add("r_error")
    return sorted(f)


# ------------------------------------------------------------------ main

def build(workers: int = 4, fresh: bool = False, out_dir: Path | None = None, subset: set[str] | None = None) -> dict:
    out_dir = out_dir or BUILD
    out_dir.mkdir(parents=True, exist_ok=True)
    ops = load_ops()
    pool = load_pool()
    if subset is not None:
        pool = [c for c in pool if c["id"] in subset]
    t0 = time.time()
    versions = r_versions()
    fp = fingerprint(versions)
    cache: dict[str, dict] = {}
    if CACHE.exists() and not fresh:
        for row in iter_jsonl(CACHE):
            if row["fp"] == fp:
                cache[row["key"]] = row["rec"]
    keyed = [(c, hashlib.sha256(f"{fp}|{case_hash(c)}".encode()).hexdigest()) for c in pool]
    todo = [c for c, k in keyed if k not in cache]
    # duplicate content hashes share one R run
    uniq: dict[str, dict] = {}
    for c, k in keyed:
        if k not in cache and k not in uniq:
            uniq[k] = c
    t_r = time.time()
    results = run_r(list(uniq.values()), ops, workers)
    r_seconds = time.time() - t_r
    new_rows = []
    for k, c in uniq.items():
        rec = dict(results[c["id"]])
        rec.pop("id", None)
        cache[k] = rec
        new_rows.append({"fp": fp, "key": k, "rec": rec})
    if new_rows:
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        with CACHE.open("w" if fresh and subset is None else "a") as fh:
            for row in new_rows:
                fh.write(json.dumps(row, separators=(",", ":")) + "\n")
    resolved, expected = [], []
    for c, k in keyed:
        rec = cache[k]
        h = case_hash(c)
        base = {"id": c["id"], "source": c["source"], "op": c["op"], "hash": h}
        if rec["status"] != "ok":
            expected.append({**base, "status": "excluded", "reason": rec["reason"], "detail": rec.get("detail", ""),
                             "warnings": rec.get("warnings", [])})
            continue
        errored = rec["result"]["kind"] == "error"
        resolved.append({**base, "inputs": rec["inputs"], "call": c["call"]})
        exp = {**base, "status": "ok", "result": rec["result"], "r_class": rec.get("r_class"),
               "domain": rec.get("domain", "other"), "warnings": rec.get("warnings", []),
               "features": features_of(c, rec["inputs"], errored)}
        if "sorted_by_query" in rec:
            exp["sorted_by_query"] = rec["sorted_by_query"]
        if "universe" in rec:
            exp["universe"] = rec["universe"]
        expected.append(exp)
    write_jsonl(out_dir / "resolved.jsonl", resolved)
    write_jsonl(out_dir / "expected.jsonl", expected)
    meta = {"versions": versions, "oracle_R_sha256": oracle_digest(), "cases": len(pool), "ok": len(resolved),
            "excluded": len(pool) - len(resolved), "r_ran": len(uniq), "r_seconds": round(r_seconds, 1),
            "wall_seconds": round(time.time() - t0, 1), "workers": workers}
    (out_dir / "oracle-meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    if subset is None and len(uniq) >= len(pool) // 2:  # a (nearly) full R run: keep its timing for the report
        RESULTS.mkdir(exist_ok=True)
        (RESULTS / "oracle-timing.json").write_text(json.dumps(
            {"versions": versions, "cases": len(pool), "r_ran": len(uniq), "r_seconds": meta["r_seconds"],
             "workers": workers}, indent=1) + "\n")
    return meta


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    ap.add_argument("--fresh", action="store_true", help="ignore the cache")
    ap.add_argument("--out", type=Path, default=None)
    ns = ap.parse_args()
    meta = build(ns.workers, ns.fresh, ns.out)
    print(json.dumps(meta, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
