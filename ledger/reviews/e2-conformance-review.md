# E2 conformance review

HEAD: `834896a`. **Changes requested: 3 major, 2 minor; no blockers.**
Paths below are relative to `experiments/e2-conformance/`.

## Findings

1. **Major — `src/ranges_e2/iranges.py:110–128`: public accessors break the promised immutability.** `x.width[0] = -5` succeeds, bypasses validation, and changes a GRanges constructed from `x`. `x.mcols.replace_column(...)` likewise changes both objects; constructors retain the caller's frame when its height matches. The immutability test only checks assignment to an attribute. **Fix:** isolate mutable constructor inputs and accessor results (array copies/read-only storage and cloned DataFrame wrappers); add mutation/aliasing tests for IRanges, GRanges, and derived objects.

2. **Major — `src/ranges_e2/_mcols.py:38–50`: concatenation silently truncates metadata.** The first frame determines each column's dtype, and later frames are cast to it. Concatenating scores `[1]` (integer) and `[2.75]` (double) produces `[1, 2]`; reversing the inputs preserves `2.75`. This violates the stated row-binding/metadata-preservation behavior. **Fix:** resolve a lossless common dtype before concatenation, or reject incompatible types rather than truncating. Add mixed-dtype metadata cases in both input orders for both range classes.

3. **Major — `src/ranges_e2/conformance/codec.py:207–218,271–273`: metadata serialization only supports a subset of its declared format.** A DFrame column `groups` containing the FORMAT.md IntegerList `[[1,2]]` decodes but cannot re-encode: `TypeError: cannot encode column groups of type List(Int32)`. The same content works only when named `revmap`. A Factor column is rejected during decode. Both are declared serialized types; the R codec accepts them as DFrame columns. Existing golden roundtrips do not cover these combinations. **Fix:** preserve metadata column type information independently of its name, separate revmap index conversion from generic list values, and add language-neutral cases plus codec roundtrips for non-revmap IntegerList and Factor columns.

4. **Minor — `conformance/r/run.R:112`: the R runner ignores the declared comparison mode.** It always calls raw `tree_equal`, while the Python runner honors `hits-query-grouped`. Thus a permitted within-query subject reorder fails only the R runner. A direct codec probe with identical Hits except `to=[2,1]` versus `[1,2]` confirms raw R comparison returns FALSE. **Fix:** resolve the case/file comparison mode in the R runner and normalize Hits accordingly; test that both runners accept within-query reorder and reject query-order changes.

5. **Minor — `conformance/cases/granges-nearest.json:7`: formal clauses contradict other checked-in clauses and behavior.** `GR.NEAREST.ALL` promises every tied candidate, but line 5 and `gr.nearest.all.tie-precede-side-only` explicitly retain only the precede side. Also, `granges-setops.json:13` says compatible pintersect pairs have `hit=TRUE` regardless of overlap, contradicting `GR.PINTERSECT.NOHIT` and its compatible/non-overlapping case. `iranges-overlaps.json:12` promises subject-sorted Hits although FORMAT.md/NOTES deliberately waive that order. **Fix:** reconcile these clause statements with the existing cases, regenerate goldens per the maintenance rule, and rebuild REPORT.md.

## Scope and approach

- All TARGET operation families have public methods; no entirely missing required operation found. Metadata preservation/serialization is partial as above.
- Suite-first approach is supported by current-branch history: format/R codec, then all 17 case groups and R goldens, precede the Python core commit. Later differential fuzzing supplements rather than replaces the suite; changes are documented.
- Four intentional xfails and known nearest/directional limits are recorded. They were not independently adjudicated here.
- Kernels are separated from objects and the harness; maintenance instructions and reports are usable. Mutable ownership and dtype handling need correction.

## Validation

From the experiment directory, with bytecode/cache writing disabled:

- `uv run pytest -p no:cacheprovider`: **1997 passed, 4 xfailed**.
- `uv run python -m ranges_e2.conformance.runner`: **901 pass, 4 xfail / 905**.
- `Rscript conformance/r/selftest.R`: **15 fixtures passed**.
- `Rscript conformance/r/run.R validate`: **905 valid cases**.
- `Rscript conformance/r/run.R check <all 17 explicit group names>`: **905 cases, 0 not passing**. Explicit groups avoid rewriting the tracked R report.
- Report builder output equals stored REPORT.md.
- Focused Python probes reproduced findings 1–3; the R codec-only probe supports finding 4. No independent R/Python fidelity hunt or fuzz run.
- Worktree unchanged; no staged files. No source/tests edited or commits created.
