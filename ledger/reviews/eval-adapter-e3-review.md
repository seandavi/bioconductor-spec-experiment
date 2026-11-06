## Findings

1. **High — NA/default conflation (thinness + translation)** — `eval/adapters/e3.py:65`. A singleton NA becomes `None`, the port's absent-argument sentinel. For example, `IRanges(start=NA, end=7, width=3)` silently derives start=5; passing a masked NA instead raises. Likewise `narrow(width=NA)` returns unchanged ranges instead of the port's NA rejection. This changes what the port computes; affects constructors, narrow, gaps and coverage. **Fix:** preserve singleton NA as a masked array; reserve `None` for R NULL. Do not implement R's NA/default behavior in the adapter.

2. **Medium — Unsupported double subscripts dispatched as supported** — `eval/adapters/e3.py:169–172`. Both integer and double vectors are shifted, although `ops.json` marks only integer subscripts index-valued and CORE.9 documents integer sequences. `i={dbl:[1.0]}` reaches the port as a float array and raises `IndexError`; the same R position expressed as integer works. **Fix:** declare double subscripts `Unsupported("subscript:double")`; do not truncate/coerce them toward R.

3. **Medium — Native atomic metadata results incorrectly unsupported** — `eval/adapters/e3.py:108–117`. `_to_column` rejects numpy arrays produced by public accessors. Thus `IRanges(start=1, width=1, w=width(x))` is unsupported, although `pl.Series("w", x.width)` supplies valid metadata to the public constructor. **Fix:** translate native ndarray/masked-array columns to dtype-preserving polars Series, preserving masks as nulls.

4. **Medium — Factor seqnames unnecessarily unsupported** — `eval/adapters/e3.py:246–247`. `build_Factor` returns an iterable string-valued `pl.Enum` Series, which public `GRanges(seqnames, ranges)` already accepts directly. The adapter blocks it before the port runs (resolved case `e2/granges-core/gr.ctor.factor-seqnames-levels`). **Fix:** remove the guard and pass the Series through; let any loss of unused factor levels remain the port's result, not an adapter repair.

5. **Medium — Public result lengths unnecessarily unsupported** — `eval/adapters/e3.py:275–276`. `length(findOverlaps(x,x))` and `length(coverage(x))` are blocked, although public `Hits.__len__` and `Rle.__len__` provide their lengths. Probes returned 1 and 7 directly versus `Unsupported` through the adapter. **Fix:** admit `Hits` and `Rle` in `op_length`.

No sorting, deduplication, clipping, overlap re-derivation, private-API use, or caught-port-error substitutions found. Static seqinfo-setter unsupported declarations match the public API. Documented integer-index/Hits shifts, strand/Seqinfo input construction and normal projection shapes checked; roundtrip passed.

## Validation

- `cd eval && just roundtrip e3`: passed; **5,866 distinct inputs**, no failures.
- `cd eval && just adapter-check e3`: passed; **9,650 cases**: 6,717 value, 1,401 port_error, 1,532 unsupported, **0 adapter_error**.
- `/tmp/e3_adapter_review_probes.py`: direct-public-API probes confirmed findings 1–5; no expected values read.
- Both worktrees remain clean; nothing staged. No source edits or commits. No `just run`, `just report`, port tests, or `build/expected.jsonl` accessed. Adapter-check generated its normal ignored raw output.

```acceptance-report
{
  "criteriaSatisfied": [
    {"id": "criterion-1", "status": "satisfied", "evidence": "Completed read-only adapter review against ADAPTERS.md and public e3 API; no source edits or commits."},
    {"id": "criterion-2", "status": "satisfied", "evidence": "Five numbered findings with locations, concrete fixes and public-API probes; both requested checks ran."}
  ],
  "changedFiles": ["/Users/davsean/.pi/agent/sessions/--Users-davsean-Documents-git-bioconductor-spec-experiment--/subagent-artifacts/outputs/f6e1b6f5-4eb8-47d4-b10c-98bc217803e6/reviews/eval-adapter-e3-review.md"],
  "testsAddedOrUpdated": [],
  "commandsRun": [
    {"command": "cd eval && just roundtrip e3", "result": "passed", "summary": "5866 distinct inputs; no failures"},
    {"command": "cd eval && just adapter-check e3", "result": "passed", "summary": "9650 cases; zero adapter_error"},
    {"command": "PYTHONDONTWRITEBYTECODE=1 uv run --frozen --project /Users/davsean/Documents/git/bioconductor-spec-experiment.worktrees/e3-spec-adr/experiments/e3-spec-adr python /tmp/e3_adapter_review_probes.py", "result": "passed", "summary": "Observed NA conflation, double-subscript rejection, native metadata rejection, factor guard and result-length guards"},
    {"command": "git status --porcelain; git diff --cached --name-only; git -C /Users/davsean/Documents/git/bioconductor-spec-experiment.worktrees/e3-spec-adr status --porcelain", "result": "passed", "summary": "Empty output: both worktrees clean; no staged files"}
  ],
  "validationOutput": ["roundtrip: 5866 distinct inputs", "9650 cases: value=6717, port_error=1401, unsupported=1532, adapter_error=0"],
  "residualRisks": ["No scored results or expected projections examined; score impact is not quantified."],
  "noStagedFiles": true,
  "diffSummary": "Review artifact only; no repository source changes.",
  "reviewFindings": ["high: e3.py:65 - singleton NA becomes absent/default", "medium: e3.py:169 - unsupported double subscripts dispatched", "medium: e3.py:108 - ndarray metadata rejected", "medium: e3.py:246 - supported factor iterable blocked", "medium: e3.py:275 - public Hits/Rle lengths blocked"],
  "manualNotes": "Recommended next step: address the findings, rerun the two adapter checks, then repeat the thinness review before scoring."
}
```
