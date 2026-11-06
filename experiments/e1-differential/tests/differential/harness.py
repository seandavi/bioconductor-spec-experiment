"""Case loading, R execution + cache, Python execution, comparison."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

import adapter
import codec

HERE = Path(__file__).parent
CASES_DIR = HERE / "cases"
CACHE_DIR = HERE / "r_cache"
R_SCRIPT = HERE / "r" / "run_cases.R"

_CASE_KEYS = ("inputs", "call", "raw")


@dataclass(frozen=True)
class Case:
    id: str
    file: str  # yaml stem
    spec: dict[str, Any]

    @property
    def r_spec(self) -> dict[str, Any]:
        """What R sees (and what the cache key hashes)."""
        return {"id": self.id, **{k: self.spec[k] for k in _CASE_KEYS if k in self.spec}}

    @property
    def hash(self) -> str:
        return hashlib.sha256(json.dumps(self.r_spec, sort_keys=True).encode()).hexdigest()


def load_cases() -> list[Case]:
    cases: list[Case] = []
    seen: dict[str, str] = {}
    for path in sorted(CASES_DIR.glob("*.yaml")):
        doc = yaml.safe_load(path.read_text()) or {}
        for spec in doc.get("cases", []):
            if spec["id"] in seen:
                raise ValueError(f"duplicate case id {spec['id']!r} in {path.name} and {seen[spec['id']]}")
            if ("call" in spec) == ("raw" in spec):
                raise ValueError(f"case {spec['id']!r}: exactly one of 'call' or 'raw'")
            seen[spec["id"]] = path.name
            cases.append(Case(spec["id"], path.stem, spec))
            if "universe" in spec:  # R select="all" result, used by compare: arbitrary
                u = {"inputs": spec.get("inputs", {}), "call": spec["universe"]}
                cases.append(Case(spec["id"] + "#universe", path.stem, u))
    return cases


# ---- R side ----------------------------------------------------------------


def r_available() -> bool:
    return shutil.which("Rscript") is not None and not os.environ.get("RANGES_E1_NO_R")


_versions: str | None = None


def r_versions() -> str:
    global _versions
    if _versions is None:
        out = subprocess.run(
            ["Rscript", "-e",
             'suppressWarnings(suppressMessages({library(IRanges);library(GenomicRanges)}));'
             'cat("R", as.character(getRversion()), "; IRanges", as.character(packageVersion("IRanges")),'
             '"; GenomicRanges", as.character(packageVersion("GenomicRanges")),'
             '"; S4Vectors", as.character(packageVersion("S4Vectors")))'],
            capture_output=True, text=True, check=True,
        ).stdout
        scripts = hashlib.sha256((R_SCRIPT.read_bytes() + (R_SCRIPT.parent / "encode.R").read_bytes())).hexdigest()[:12]
        _versions = f"{out.strip()}; harness {scripts}"
    return _versions


def _cache_path(file: str) -> Path:
    return CACHE_DIR / f"{file}.json"


def _load_cache(file: str) -> dict[str, Any]:
    p = _cache_path(file)
    return json.loads(p.read_text()) if p.exists() else {"versions": None, "entries": {}}


_fresh: dict[str, dict[str, Any]] = {}  # file -> cache doc validated this session


def refresh_cache(file: str, cases: list[Case], *, force: bool) -> dict[str, Any]:
    """Return the cache doc for `file`, re-running R for stale entries if possible."""
    if file in _fresh:
        return _fresh[file]
    doc = _load_cache(file)
    ver_ok = True
    if r_available():
        ver_ok = doc.get("versions") == r_versions()
    stale = [
        c for c in cases
        if force or not ver_ok or doc["entries"].get(c.id, {}).get("hash") != c.hash
    ]
    if stale and r_available():
        results, ver = run_r(stale)
        if not ver_ok or force:
            doc = {"versions": ver, "entries": {}}
        doc["versions"] = ver
        for c in stale:
            doc["entries"][c.id] = {"hash": c.hash, "result": results[c.id]}
        live = {c.id for c in cases}
        doc["entries"] = {k: doc["entries"][k] for k in sorted(doc["entries"]) if k in live}
        CACHE_DIR.mkdir(exist_ok=True)
        _cache_path(file).write_text(json.dumps(doc, indent=1) + "\n")
    _fresh[file] = doc
    return doc


def run_r(cases: list[Case]) -> tuple[dict[str, Any], str]:
    with tempfile.TemporaryDirectory() as td:
        tdp = Path(td)
        (tdp / "in.json").write_text(json.dumps([c.r_spec for c in cases]))
        proc = subprocess.run(
            ["Rscript", str(R_SCRIPT), str(tdp / "in.json"), str(tdp / "out.json")],
            capture_output=True, text=True,
        )
        if proc.returncode != 0 or not (tdp / "out.json").exists():
            raise RuntimeError(f"Rscript failed:\n{proc.stderr[-2000:]}")
        out = json.loads((tdp / "out.json").read_text())
    return out["results"], out["versions"]


def r_result(case: Case, all_cases: list[Case], *, force: bool = False) -> dict[str, Any] | None:
    """Cached (or freshly computed) R result; None if stale and R unavailable."""
    same_file = [c for c in all_cases if c.file == case.file]
    doc = refresh_cache(case.file, same_file, force=force)
    entry = doc["entries"].get(case.id)
    if entry is None or entry["hash"] != case.hash:
        return None
    return entry["result"]


# ---- Python side -----------------------------------------------------------


_BUGS = (TypeError, AttributeError, NameError, AssertionError, NotImplementedError, ZeroDivisionError)


def py_result(case: Case) -> dict[str, Any]:
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            value = adapter.encode(adapter.run_python(case.spec))
            out: dict[str, Any] = {"value": value}
        except _BUGS as e:  # a Python defect, never a legitimate R-style error
            out = {"error": f"{type(e).__name__}: {e}", "bug": True}
        except Exception as e:  # noqa: BLE001 - ValueError/IndexError/... = R-style error
            out = {"error": f"{type(e).__name__}: {e}"}
    out["warnings"] = [str(x.message) for x in w]
    return out


# ---- comparison ------------------------------------------------------------


def _arbitrary_ok(r_all: dict[str, Any], py: dict[str, Any]) -> str | None:
    """Python select='arbitrary' picks must be valid hits of R's select='all'."""
    hits = r_all["value"]
    per_q: dict[int, set[int]] = {}
    for f, t in zip(hits["from"], hits["to"]):
        per_q.setdefault(f, set()).add(t)
    vals = py["value"]["values"]
    for q, pick in enumerate(vals, start=1):
        valid = per_q.get(q)
        if pick is None and valid:
            return f"query {q}: NA but R has hits {sorted(valid)}"
        if pick is not None and (not valid or pick not in valid):
            return f"query {q}: picked {pick}, R hits {sorted(valid or [])}"
    return None


def evaluate(case: Case, r: dict[str, Any], py: dict[str, Any], universe: dict[str, Any] | None = None) -> str | None:
    """None if the case passes, else a failure description."""
    if py.get("bug"):
        return f"Python raised an internal error (not an R-style error): {py['error']}"
    div = case.spec.get("divergence")
    if div is not None:
        want = div["python"]
        got = {"error": True} if "error" in py else py["value"]
        d = codec.first_diff(want, got)
        if d:
            return f"python no longer matches ADR-{div['adr']} value: {d}"
        return None
    if case.spec.get("compare") == "arbitrary":
        if "error" in py or "error" in r:
            return None if ("error" in py and "error" in r) else "error on one side only"
        assert universe is not None
        return _arbitrary_ok(universe, py)
    r_err, p_err = "error" in r, "error" in py
    if r_err or p_err:
        if r_err and p_err:
            return None
        return f"R {'errors: ' + r['error'] if r_err else 'ok'}; Python {'errors: ' + py['error'] if p_err else 'ok'}"
    fn = case.spec.get("call", {}).get("fn")
    d = codec.first_diff(adapter.normalize_r(fn, r["value"]), py["value"])
    if d:
        return d
    if case.spec.get("check_warnings") and bool(r["warnings"]) != bool(py["warnings"]):
        return f"warnings differ: R {r['warnings']} vs Python {py['warnings']}"
    return None
