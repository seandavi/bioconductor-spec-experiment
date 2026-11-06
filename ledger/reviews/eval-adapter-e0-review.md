1. **High — thinness / silent deduplication:** `eval/adapters/e0.py:275,477`. Both `dict(zip(...))` conversions silently discard repeated names, keeping the last value. Probe: `seqlengths=c(a=5,a=8)` becomes `{'a': 8}`; `coverage(width=c(a=3,a=4))` becomes `{'a': 4}`. **Fix:** reject duplicate names as `Unsupported` before either mapping conversion, as `mcols_dict` already does. The public mapping API cannot represent these vectors losslessly.

2. **High — thinness / NA becomes a default:** `eval/adapters/e0.py:75-79,270`. Scalar character NA becomes Python `None`; `GRanges(..., strand=None)` computes `'*'`, whereas its public vector form `strand=[None]` raises `ValueError`. R's public constructor converts missing strand to `'*'` with a warning, so this adapter conversion masks a real port deviation. **Fix:** translate atomic `strand` arguments as lists, including length-one vectors; reserve Python `None` for R NULL. Do not repair the port's NA behavior.

3. **Medium — named-width translation depends on length:** `eval/adapters/e0.py:473`. Only multi-element named widths become mappings; a named singleton loses its name through `to_arg`. Probe with seqlevels `a,b`: `width=c(a=3)` produces lengths `{a:3,b:3}`, while the public mapping form `width={'a':3}` produces `{a:3,b:1}`. **Fix:** preserve named-width mappings for singleton vectors too; apply scalar conversion to unnamed scalar widths. R's public API rejects this incomplete named width rather than treating it as an unnamed scalar. Preserve the port's mapping behavior; do not add R validation.

4. **Medium — supported Factor metadata declared unsupported:** `eval/adapters/e0.py:107-109`. `build_Factor` produces a `pl.Series` with Enum dtype, but constructors reject that same native value in `...` as `mcols: non-atomic column`. Direct public `IRanges(..., mcols={'f': factor_series})` succeeds and retains Enum dtype. GRanges uses the same rejecting helper. **Fix:** accept supported native `pl.Series` columns alongside RVector columns, preserving dtype and the requested column name.

Checks:
- `cd eval && just roundtrip e0`: **failed**, 5,866 distinct inputs; four unsupported NA-name inputs (three IRanges, one GRanges), no reported projection mismatches.
- `cd eval && just adapter-check e0`: **passed**, 9,650 cases; 7,650 values, 1,323 port errors, 677 unsupported, zero adapter errors. This does not establish semantic correctness.
- `/tmp/e0-review-probes.py`: reproduced all four findings using the port's public API. Ad-hoc `Rscript -e` public-API probes confirmed NA-strand conversion and named-width rejection; no oracle/expected files used. No port tests, goldens, scored results, or `build/expected.jsonl` read; no `just run`/`just report` invoked.
- No worktree source edits or commits; final git status clean and staged-file list empty. Static unsupported setter declarations match the absence of public setters. No additional sorting, clipping, result re-derivation, or caught-port-error substitution found.

```acceptance-report
{
  "criteriaSatisfied": [
    {"id":"criterion-1","status":"satisfied","evidence":"Completed read-only adapter review; only requested external review artifact written."},
    {"id":"criterion-2","status":"satisfied","evidence":"Four numbered findings with locations, concrete fixes, public-API probes, and both requested command results."}
  ],
  "changedFiles": [],
  "testsAddedOrUpdated": [],
  "commandsRun": [
    {"command":"cd eval && just roundtrip e0","result":"failed","summary":"5866 inputs; four unsupported NA-name inputs."},
    {"command":"cd eval && just adapter-check e0","result":"passed","summary":"9650 cases; zero adapter_error; 677 unsupported."},
    {"command":"PYTHONDONTWRITEBYTECODE=1 <e0-control>/.venv/bin/python /tmp/e0-review-probes.py","result":"passed","summary":"Reproduced duplicate-name loss, scalar NA/default collision, singleton width name loss, and rejected supported Factor metadata."},
    {"command":"Rscript -e (ad-hoc GenomicRanges public-API probes)","result":"passed","summary":"R converts NA strand to wildcard with warning; rejects named width c(a=3L) for seqlevels a,b."},
    {"command":"git status --short; git diff --cached --name-only","result":"passed","summary":"Both outputs empty."}
  ],
  "validationOutput": ["See numbered findings and command summaries above."],
  "residualRisks": ["Findings remain unfixed; review-only assignment.","Ad-hoc R checks used locally installed packages, which warn they were built under R 4.6.1; no full scoring performed."],
  "noStagedFiles": true,
  "diffSummary": "No source diff; external review artifact only.",
  "reviewFindings": ["High: duplicate named-vector entries silently discarded.","High: explicit scalar NA strand becomes wildcard default.","Medium: singleton named coverage width loses its name.","Medium: supported Factor metadata rejected."],
  "manualNotes": "Recommended next step: fix findings and rerun the two adapter checks; independent reviewer acceptance still required."
}
```
