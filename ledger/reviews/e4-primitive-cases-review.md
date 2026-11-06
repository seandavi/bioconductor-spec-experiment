# E4 review — HEAD `59de392`

## Scope and validation

- All TARGET operation families have public APIs and shared cases. No family is missing. Nearest/precede/follow require an explicit subject in both classes; README documents this limit only for GRanges.
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -p no:cacheprovider`: **2250 passed, 6 skipped**. Cache writes disabled for read-only review.
- `Rscript runners/r/run_cases.R --check`: **goldens up to date**. Package-build-version warnings only.
- Approach followed: primitive YAML constructor inputs, clause IDs, native R execution, generated projected JSON, Python projections, ADR-linked exceptions, staged commits. Comparator defects below weaken conformance guarantees.
- Design: small kernels and grouping helpers are readable. Boundary validation and result-kind inference need strengthening.
- No independent R/Python fidelity search. Probes test this experiment's own contracts only. Worktree and index remain clean.

All paths below are relative to `experiments/e4-primitive-cases/`.

## Findings

1. **Major — “exact” comparison accepts unequal coordinates.** `tests/conformance/compare.py:32–33`
   `_eq` applies numeric tolerance even to integers. Comparing vector `[2147483647]` against `[2147483646]` returns no problems. Thus an off-by-one coordinate can pass without an ADR, contrary to PLAN §4.3's exact mode. Unequal floats also pass exact comparison.
   **Fix:** use equality for exact mode, always exact for integer coordinates/indices. If float tolerance is needed, introduce an explicit comparison mode and ADR. Add comparator regression tests for both mismatches.

2. **Major — result construction incompletely enforces int32 bounds.** `src/ranges_e4/iranges.py:93–98`
   `_init` checks only maximum end and minimum start. `IRanges(start=[2147483647], width=[0]).shift(1)` succeeds with start `2147483648`, violating ADR-0006 and IR-EDGE-2. `IRanges(start=[-2147483647], width=[0])` accepts end `-2147483648`. Inferred width can also exceed int32: start `-2147483647`, end `2147483647` yields width `4294967295` despite PLAN's integer-storage contract.
   **Fix:** centrally validate both bounds of start/end and the width limit, including inferred coordinates and `_make` results. Add constructor and operation boundary tests for zero-width ranges in both classes.

3. **Major — GRanges coverage bypasses vector-length validation.** `src/ranges_e4/granges.py:498–515`
   Weight is indexed by group before its global length is checked. For two ranges on `a,b`, `coverage(weight=[1,2,3])` silently uses `[1,2]`; `coverage(weight=[])` raises `IndexError`. ADR-0007 requires non-scalar length mismatches to raise `ValueError`.
   **Fix:** validate numeric weight and its length against the whole input before splitting by sequence. Reuse this validation for circular coverage. Add too-short, too-long, empty, and valid-vector tests.

4. **Major — distance rejects documented zero-length inputs.** `src/ranges_e4/iranges.py:320–324`; `src/ranges_e4/granges.py:473–478`
   IR-DISTANCE-3 and GR-DISTANCE-4 promise zero-length output for zero-length inputs. Both classes instead raise `ValueError` for lengths `(0,2)` and `(2,0)` because length checking precedes empty handling. Existing cases exercise only `(0,0)`.
   **Fix:** return an empty masked integer result when either input is empty, before recycling checks or broadcasting. Add asymmetric empty-input cases and unit tests.

5. **Minor — advertised non-exact modes can hide unrelated errors.** `tests/conformance/compare.py:37–38,53–55`
   `set` recursively sorts the outer list, so revmaps `[[1],[2]]` and `[[2],[1]]` compare equal although PLAN §4.3 permits ignoring only member order per element. `multiset` uses truncating zip: start `[1,99]`, width `[1]` compares equal to start `[1]`, width `[1]`. Neither mode is currently exercised by shared cases.
   **Fix:** preserve outer row order for list-field set comparison; validate row-field lengths before multiset comparison. Add direct positive and negative comparator tests.

6. **Minor — circular-coverage clauses contradict each other.** `clauses.md:307`
   GR-EDGE-1 says coverage does not wrap circular sequences. GR-EDGE-6, circular coverage cases, implementation, and `test_circular_coverage_wraps` say it does.
   **Fix:** retire/correct GR-EDGE-1's coverage claim under the document's append-only convention; retain its intra-range/reduce assertions and reference GR-EDGE-6 for coverage.

7. **Minor — empty coverage loses its projected result kind.** `tests/conformance/project.py:72–78`
   `project(GRanges([], start=[], width=[]).coverage())` produces `kind: vector`, not the documented `RleList`: an empty coverage dict fails the truthiness guard and becomes a named-vector dict. Current empty coverage cases retain seqlevels, hiding this ambiguity.
   **Fix:** pass expected result-kind context from the operation dispatcher for coverage, including empty dicts. Add a shared coverage case with no ranges and no seqlevels; ensure empty named-vector accessors still project as vectors.

## Reproduction

From the experiment directory:

```sh
PYTHONDONTWRITEBYTECODE=1 uv run python /tmp/e4_primitive_cases_review_01fcc0b8.py
```

The permitted `/tmp` probe reproduces findings 1–5 and 7 without writing repository files. Finding 6 is a direct artifact contradiction. No source changes or tests were added during review.
