"""`just report`: matrix port x family x source x outcome, per-source headline, mismatch kinds, features, exclusions."""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict

from .common import BUILD, CASES, RESULTS, iter_jsonl, load_families, ports, read_jsonl
from .compare import OUTCOMES

SOURCE_ORDER = ["heldout", "e1", "e1-redteam", "e2", "e4"]
HOME = {"e0": set(), "e0b": set(), "e1": {"e1", "e1-redteam"}, "e2": {"e2"}, "e3": set(), "e4": {"e4"}}
SHORT = {"pass": "pass", "mismatch": "mism", "adapter_error": "adapt", "unsupported": "unsup",
         "port_error_unexpected": "perr"}


def pct(n: int, d: int) -> str:
    return "-" if d == 0 else f"{100 * n / d:.1f}%"


def md_table(header: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def counts_cell(c: Counter) -> str:
    total = sum(c.values())
    if total == 0:
        return "-"
    extra = " ".join(f"{SHORT[o]}{c[o]}" for o in OUTCOMES[1:] if c[o])
    return f"{c['pass']}/{total} ({pct(c['pass'], total)})" + (f" [{extra}]" if extra else "")


def build(results=RESULTS, build_dir=BUILD, cases=CASES, registry=None, selftest_ports=None) -> tuple[str, dict]:
    """Render the report from scored results; paths and registry are parameters so tests can use a slice."""
    families = load_families()
    registry = registry or ports()
    expected = {r["id"]: r for r in iter_jsonl(build_dir / "expected.jsonl")}
    oracle_meta = json.loads((build_dir / "oracle-meta.json").read_text())
    unconv = read_jsonl(cases / "unconvertible.jsonl")
    conv = json.loads((cases / "conversion-summary.json").read_text())

    scored: dict[str, list[dict]] = {}
    meta: dict[str, dict] = {}
    for name in registry:
        f = results / "scored" / f"{name}.jsonl"
        if f.exists():
            scored[name] = read_jsonl(f)
            m = results / "scored" / f"{name}.meta.json"
            meta[name] = json.loads(m.read_text()) if m.exists() else {}
    real = [n for n in scored if not registry[n].get("selftest")]
    selftest = [n for n in scored if registry[n].get("selftest")]

    # matrix: port -> family -> source -> outcome -> n
    matrix: dict = {}
    for name, rows in scored.items():
        mx: dict = defaultdict(lambda: defaultdict(Counter))
        for r in rows:
            mx[families[r["op"]]][r["source"]][r["outcome"]] += 1
        matrix[name] = {f: {s: dict(c) for s, c in d.items()} for f, d in mx.items()}

    L: list[str] = ["# Fidelity to R: evaluation report", ""]
    timing_file = results / "oracle-timing.json"
    timing = json.loads(timing_file.read_text()) if timing_file.exists() else None
    L += [f"R: `{oracle_meta['versions']}`", f"Oracle: {oracle_meta['cases']} cases, {oracle_meta['ok']} scored, "
          f"{oracle_meta['excluded']} excluded."
          + (f" Full R run: {timing['r_ran']} distinct cases in {timing['r_seconds']} s ({timing['workers']} workers)."
             if timing else ""), ""]

    # --- pool reconciliation
    L += ["## Pool and reconciliation", ""]
    exc = Counter()
    exc_src = Counter()
    for e in expected.values():
        if e["status"] != "ok":
            exc[e["reason"]] += 1
            exc_src[e["source"]] += 1
    ok_src = Counter(e["source"] for e in expected.values() if e["status"] == "ok")
    rows = []
    for suite in ("e1", "e2", "e4"):
        c = conv[suite]
        srcs = [s for s in SOURCE_ORDER if s.startswith(suite)]
        n_ok = sum(ok_src[s] for s in srcs)
        n_ex = sum(exc_src[s] for s in srcs)
        assert c["read"] == c["converted"] + c["unconvertible"] and c["converted"] == n_ok + n_ex, suite
        rows.append([suite, c["read"], c["converted"], c["unconvertible"], n_ex, n_ok])
    h_ok, h_ex = ok_src["heldout"], exc_src["heldout"]
    rows.append(["heldout", h_ok + h_ex, h_ok + h_ex, 0, h_ex, h_ok])
    L += [md_table(["suite", "read/generated", "converted", "unconvertible", "R-excluded", "scored"], rows), "",
          "read = converted + unconvertible; converted = R-excluded + scored (asserted when this report is built).", ""]
    if unconv:
        L += ["Unconvertible (listed in `cases/unconvertible.jsonl`):", ""]
        uc = Counter((u["source"], u["reason"]) for u in unconv)
        L += [md_table(["suite", "reason", "n"], [[s, r, n] for (s, r), n in sorted(uc.items())]), ""]
    if exc:
        L += ["R-excluded cases: " + ", ".join(f"{r} x{n}" for r, n in exc.most_common()) + ".", ""]

    # --- headline
    if real:
        L += ["## Headline: pass rate per port x source", "",
              "Pass rate = pass / scored; unsupported, adapter errors and port errors stay in the denominator. "
              "**heldout is the headline** (no port has seen it). `*` = the port's home suite (it was built against it).", ""]
        header = ["port"] + SOURCE_ORDER + ["all"]
        rows = []
        for name in real:
            by = defaultdict(Counter)
            for r in scored[name]:
                by[r["source"]][r["outcome"]] += 1
            row = [name]
            for s in SOURCE_ORDER:
                c = by[s]
                row.append(("*" if s in HOME.get(name, set()) else "") + (counts_cell(c) if c else "-"))
            tot = Counter()
            for c in by.values():
                tot.update(c)
            row.append(f"{tot['pass']}/{sum(tot.values())} ({pct(tot['pass'], sum(tot.values()))})")
            rows.append(row)
        L += [md_table(header, rows), "",
              "Cell: `pass/scored (rate) [mism = mismatch, adapt = adapter_error, unsup = unsupported, perr = port_error_unexpected]`.", ""]

        # --- matrices per port
        for name in real:
            L += [f"## {name}: family x source", ""]
            fams = list(dict.fromkeys(families.values()))
            rows = []
            for f in fams:
                row = [f]
                for s in SOURCE_ORDER:
                    row.append(counts_cell(Counter(matrix[name].get(f, {}).get(s, {}))))
                rows.append(row)
            L += [md_table(["family"] + SOURCE_ORDER, rows), ""]

            # mismatch kinds
            kinds: dict = defaultdict(Counter)
            rules = Counter()
            for r in scored[name]:
                if r["outcome"] == "mismatch":
                    for k in r["detail"]["kinds"]:
                        kinds[r["source"]][k] += 1
                    rules[r["rule"]] += 1
            if kinds:
                ks = sorted({k for c in kinds.values() for k in c})
                L += [f"Mismatch kinds (a case counts once per kind): rule used {dict(rules)}", "",
                      md_table(["source"] + ks, [[s] + [kinds[s].get(k, 0) for k in ks] for s in SOURCE_ORDER if s in kinds]), ""]

            # features (heldout)
            feat: dict = defaultdict(Counter)
            for r in scored[name]:
                if r["source"] == "heldout":
                    for f in expected[r["id"]]["features"]:
                        feat[f][r["outcome"]] += 1
            if feat:
                L += [f"{name} heldout by input feature:", "",
                      md_table(["feature", "pass/scored"], [[f, counts_cell(feat[f])] for f in sorted(feat)]), ""]

        # --- runs
        L += ["## Runs", "",
              md_table(["port", "tag", "frozen sha (before/after)", "clean", "adapter sha256", "run s"],
                       [[n, meta[n].get("tag"), f"{meta[n].get('sha_before')}/{meta[n].get('sha_after')}",
                         f"{meta[n].get('clean_before')}/{meta[n].get('clean_after')}",
                         (meta[n].get("adapter_sha256") or "")[:12], meta[n].get("run_seconds")] for n in real]), ""]

    # --- self-tests
    if selftest:
        L += ["## Self-tests", ""]
        rows = []
        for name in selftest:
            tot = Counter(r["outcome"] for r in scored[name])
            rows.append([name, sum(tot.values())] + [tot[o] for o in OUTCOMES])
        L += [md_table(["adapter", "scored"] + list(OUTCOMES), rows), ""]
    return "\n".join(L) + "\n", {"oracle": oracle_meta, "matrix": matrix}


def main() -> int:
    text, data = build()
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "report.md").write_text(text)
    (RESULTS / "matrix.json").write_text(json.dumps(data, indent=1, sort_keys=True) + "\n")
    print(f"wrote {RESULTS / 'report.md'} ({len(text.splitlines())} lines), matrix.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
