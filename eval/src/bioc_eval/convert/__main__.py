"""`python -m bioc_eval.convert`: write cases/{e1,e2,e4}.jsonl, unconvertible.jsonl, conversion-summary.json."""

from __future__ import annotations

import json
import sys
from collections import Counter

from ..common import CASES, write_jsonl
from . import e1, e2, e4


def main() -> int:
    unconv_all: list[dict] = []
    summary: dict = {}
    seen: set[str] = set()
    for name, mod in (("e1", e1), ("e2", e2), ("e4", e4)):
        cases, unconv, stats = mod.convert()
        for c in cases:
            if c["id"] in seen:
                raise SystemExit(f"duplicate case id {c['id']}")
            seen.add(c["id"])
        cases.sort(key=lambda c: c["id"])
        write_jsonl(CASES / f"{name}.jsonl", cases)
        unconv_all += unconv
        by_source = Counter(c["source"] for c in cases)
        summary[name] = {
            "read": stats["read"], "converted": len(cases), "unconvertible": len(unconv),
            "converted_by_source": dict(by_source),
            "unconvertible_by_reason": dict(Counter(u["reason"] for u in unconv).most_common()),
        }
        assert stats["read"] == len(cases) + len(unconv), f"{name}: read != converted + unconvertible"
        print(f"{name}: read {stats['read']}, converted {len(cases)} {dict(by_source)}, unconvertible {len(unconv)}")
        for reason, n in Counter(u["reason"] for u in unconv).most_common():
            print(f"    {n:5d}  {reason}")
    unconv_all.sort(key=lambda u: (u["source"], u["origin"]["file"], u["origin"]["id"]))
    write_jsonl(CASES / "unconvertible.jsonl", unconv_all)
    (CASES / "conversion-summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
