"""Audit experiment child transcripts for isolation violations.

Usage: uv run --no-project python -I harness/audit_isolation.py EXP [--no-r] TRANSCRIPT...

Scans every tool-call `command`/`path` argument in pi subagent transcripts
(.jsonl) and reports:
- paths under /Users/... outside the experiment's own worktree (excluding /tmp),
- reads of background/, harness/, other worktrees or other refs,
- BiocPy packages (install, import or source reads),
- with --no-r: any R / Rscript invocation or R library read.
Exit code 1 if anything was found.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ARG_RE = re.compile(r'"(?:command|path)":\s*"((?:[^"\\]|\\.)*)"')
HOME_PATH_RE = re.compile(r"/Users/[^\s\"'\\;|&)]+")
# Module names are lowercase; case-sensitive so `import IRanges` (R-style class names) doesn't match.
BIOCPY_RE = re.compile(
    r"(pip|uv)\s+(install|add)[^\n]*\b(iranges|genomicranges|biocframe|biocutils|dolomite)"
    r"|site-packages/(iranges|genomicranges|biocframe|biocutils|dolomite)"
    r"|(^|\n)\s*(import|from)\s+(iranges|genomicranges|biocframe|biocutils)\b"
)
REFS_RE = re.compile(r"git\s+(log|show|diff|checkout|switch|branch)[^\n]*(--all|exp/e|\bmain\b|ledger)")
R_RE = re.compile(
    r"\bRscript\b|(^|[;&|(]\s*)R\s+(-e|-f|-q|--|CMD)"
    r"|/library/(IRanges|GenomicRanges|S4Vectors|BiocGenerics)/|\.Rd\b"
)


def tool_args(transcript: Path) -> list[str]:
    out: list[str] = []
    for line in transcript.read_text(errors="replace").splitlines():
        try:
            obj = json.loads(line)
        except json.JSONDecodeError:
            continue
        message = obj.get("message")
        if not isinstance(message, dict) or message.get("role") != "assistant":
            continue  # tool calls live in assistant turns; prompts mention forbidden things
        for m in ARG_RE.finditer(json.dumps(obj)):
            out.append(json.loads('"' + m.group(1) + '"'))
    return out


def audit(exp: str, transcripts: list[Path], no_r: bool) -> list[str]:
    own = f"worktrees/{exp}"
    problems: list[str] = []
    for t in transcripts:
        for arg in tool_args(t):
            short = arg[:140].replace("\n", " ")
            for p in HOME_PATH_RE.findall(arg):
                if own not in p and "/.pi/" not in p and "/.cache/uv" not in p:
                    problems.append(f"{t.name}: outside path {p}")
            if re.search(r"(^|/)background/|/harness/|worktrees/(?!" + re.escape(exp) + r")", arg):
                problems.append(f"{t.name}: forbidden dir: {short}")
            if BIOCPY_RE.search(arg):
                problems.append(f"{t.name}: BiocPy: {short}")
            if REFS_RE.search(arg):
                problems.append(f"{t.name}: other refs: {short}")
            if no_r and R_RE.search(arg):
                problems.append(f"{t.name}: R use: {short}")
    return problems


def main() -> int:
    args = sys.argv[1:]
    no_r = "--no-r" in args
    args = [a for a in args if a != "--no-r"]
    exp, files = args[0], [Path(a) for a in args[1:]]
    problems = audit(exp, files, no_r)
    for p in dict.fromkeys(problems):
        print(p)
    print(f"{exp}: {len(files)} transcripts, {len(set(problems))} findings")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
