"""Fail if any heldout seed value appears in an experiment worktree (S4 isolation check)."""

from __future__ import annotations

import os
import re
import sys

from .common import EVAL, WORKTREES, ports

SKIP_DIRS = {".git", ".venv", "__pycache__", "node_modules"}
TOKEN = re.compile(rb"\b(?:\d{17,20}|[0-9a-f]{15,16})\b")


def main() -> int:
    values = {int(s) for s in (EVAL / "heldout" / "seeds.txt").read_text().split()}
    dec = {str(v).encode() for v in values}
    hexed = {f"{v:x}".encode() for v in values}
    dirs = sorted({WORKTREES / spec["project"].split("/")[0] for spec in ports().values()
                   if not spec.get("selftest") and not spec["project"].startswith("eval")})
    hits = nfiles = 0
    for d in dirs:
        for root, subdirs, files in os.walk(d):
            subdirs[:] = [s for s in subdirs if s not in SKIP_DIRS]
            for f in files:
                path = os.path.join(root, f)
                try:
                    data = open(path, "rb").read()
                except OSError:
                    continue
                nfiles += 1
                for tok in set(TOKEN.findall(data)):
                    if tok in dec or tok in hexed:
                        hits += 1
                        print(f"SEED HIT {tok.decode()} in {path}")
    print(f"seed scan: {len(values)} seeds, {len(dirs)} worktrees, {nfiles} files, {hits} hits")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
