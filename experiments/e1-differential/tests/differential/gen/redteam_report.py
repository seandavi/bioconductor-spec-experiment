"""Render per-case findings from fresh/cached R results and current Python results."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parents[1]))
import harness  # noqa: E402

ROOT = Path(__file__).parents[3]
OUT = ROOT / "redteam"


def vec(values: list[Any]) -> str:
    return "[" + ",".join("NA" if v is None else str(v).lower() if isinstance(v, bool) else str(v) for v in values) + "]"


def summary(result: dict[str, Any]) -> str:
    if "error" in result:
        text = str(result["error"]).splitlines()[0]
        return "error: " + (text[:92] + "…" if len(text) > 92 else text)
    value = result["value"]
    if value is None:
        return "NULL"
    cls = value.get("class")
    if cls == "vector":
        return vec(value["values"])
    if cls == "Hits":
        return ", ".join(f"{q}→{s}" for q, s in zip(value["from"], value["to"])) or "no hits"
    if cls in ("IRanges", "GRanges"):
        coords = ",".join(f"{s}/{w}" for s, w in zip(value["start"], value["width"])) or "empty"
        parts = [coords]
        if value["names"] is not None:
            parts.append("names=" + vec(value["names"]))
        if value["mcols"]:
            for name, column in value["mcols"].items():
                if column["class"] == "vector":
                    parts.append(name + "=" + vec(column["values"]))
                else:
                    parts.append(name + "=" + str(column["elements"]))
        if cls == "GRanges":
            if value["start"]:
                parts.append("strand=" + vec(value["strand"]))
            if not value["start"] or len(value["seqinfo"]["seqnames"]) > 1:
                parts.append("levels=" + vec(value["seqinfo"]["seqnames"]))
        return "; ".join(parts)
    if cls == "Seqinfo":
        return "levels=" + vec(value["seqnames"]) + "; circular=" + vec(value["is_circular"])
    if cls == "Rle":
        return "runs=" + vec(value["values"]) + "; lengths=" + vec(value["lengths"])
    if cls == "RleList":
        return "; ".join(name + ": " + summary({"value": item}) for name, item in zip(value["names"], value["elements"]))
    return str(value)


def main() -> None:
    all_cases = harness.load_cases()
    OUT.mkdir(exist_ok=True)
    for group in ("circular", "metadata", "edges"):
        cases = [c for c in all_cases if c.file == f"redteam_{group}"]
        rows: list[str] = []
        failures = 0
        for case in cases:
            r = harness.r_result(case, all_cases)
            if r is None:
                raise RuntimeError(f"{case.id}: missing R result")
            py = harness.py_result(case)
            failure = harness.evaluate(case, r, py)
            if bool(failure) != bool(case.spec.get("xfail")):
                raise RuntimeError(f"{case.id}: stale xfail status: {failure}")
            failures += bool(failure)
            div = case.spec.get("divergence")
            verdict = (f"Accepted divergence (ADR-{div['adr']})" if div else
                       "ADR judgment" if case.spec.get("xfail", "").startswith("Judgment") else "Bug") if (failure or div) else "Agrees"
            p = summary(py) if (failure or div) else "same" if "value" in r else "also errors"
            rows.append(f"| `{case.id}` | {summary(r)} | {p} | {verdict} |".replace("\n", " "))
        text = f"# Red team: {group}\n\n"
        text += f"{len(cases)} cases: {len(cases) - failures} agree, {failures} strict expected failures. Accepted divergences are pinned to an ADR.\n\n"
        text += "Legend: coordinates are start/width; hits/indices are R's 1-based values; NA is missing.\n"
        text += "Tables abbreviate outputs; the comparator checks complete metadata, strand, and seqinfo.\n"
        text += "'Agrees' includes matching errors. 'ADR judgment' is unresolved, not permission to diverge.\n\n"
        text += "| Case id | R | Python | View |\n|---|---|---|---|\n" + "\n".join(rows) + "\n"
        (OUT / f"{group}.md").write_text(text)
        print(f"{group}: {len(cases)} cases, {failures} expected failures")


if __name__ == "__main__":
    main()
