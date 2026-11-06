# E3 review — `120eab8`

**Changes requested:** one major finding, two minor findings. Paths below are relative to `experiments/e3-spec-adr/`.

## Findings

1. **Major — `src/ranges_e3/_subscript.py:21–38`: missing subscripts silently select data.** CORE.9 explicitly refuses NA positions. Masked arrays enter the ordinary ndarray path; `_positions` loses their mask. For example:
   ```python
   x = IRanges(start=[10, 20], width=2)
   x[np.ma.array([1], mask=[True])].start.tolist()  # [20], should raise IndexError
   ```
   More dangerously, indexing a subject with a missing `find_overlaps(select="first")` result selects its first row. Both IRanges and GRanges are affected. Masked boolean subscripts also succeed. An empty boolean array bypasses the length check and returns empty even for nonempty input, another CORE.9 violation.
   **Fix:** reject masked entries before resolving subscripts; check boolean dimensionality/length before the generic empty-array return. Add CORE.9/GR.4 regression tests for missing integer/boolean subscripts, overlap-result subscripts, and empty boolean masks. Keep empty integer subsetting valid.

2. **Minor — `spec/90-seqinfo.md:16`: Seqinfo preservation contradicts implementation and tests.** SEQI.5 promises unchanged Seqinfo unless another clause states an exception; only concat is named. GR.12/GR.13 do not specify merging. Yet `x.union(y)`, `intersect`, `setdiff`, `pintersect`, and `psetdiff` all return merged levels: with x on `a` and y on `b`, the result levels are `["a", "b"]`, even for empty intersection. `tests/test_granges_inter.py:95–101` explicitly requires union's conflicting behavior.
   **Fix:** state the two-operand Seqinfo merge policy in GR.12/GR.13 and reference these exceptions from SEQI.5. Cover result metadata for all six set operations, including unused levels and attribute filling. This restores the spec's authority rather than letting tests define undocumented semantics.

3. **Minor — `spec/95-granges.md:47`: GR.13 gives contradictory strand rules.** Its first rule makes `psetdiff` take y's strand whenever x is `*`, including non-comparable pairs. Its later non-comparable rule says both parallel difference/intersection keep x's strand. The example on line 48 and test at `tests/test_granges_inter.py:194` implement the first rule: `a:1-8:*` minus `b:5-9:+` becomes `a:1-8:+`.
   **Fix:** say non-comparable `psetdiff` keeps x's coordinates, with strand determined by the earlier rule; only `pintersect` unconditionally keeps x's strand there.

## Scope, approach, quality

- All required operation names are public on both classes; concat and required accessors are present. No missing operations found.
- Recorded partial semantics: circular overlaps/coverage remain linear (ADR-0014); zero-width self-nearest and GRanges directional neighbors involving `[1,0]` are implementation-defined. GRanges nearest-all ties deliberately differ (ADR-0011). These are documented choices, not newly investigated R discrepancies.
- Spec-first ordering is supported by current-branch history: spec `0aa3554`, ADRs `e73d4c8`, then implementation `5af5431` onward. Staged commits carry the required attribution. Later amendments are documented. Findings 2–3 are local spec drift, not absence of the approach.
- Shared numerical kernels and GRanges grouping provide a coherent seam. Clause-citing tests and maintenance guidance are useful; citation coverage alone does not ensure semantic consistency.

## Validation

- `cd experiments/e3-spec-adr && PYTHONDONTWRITEBYTECODE=1 PYTEST_ADDOPTS='-p no:cacheprovider' uv run --frozen pytest`: **162 passed**.
- Read-only Python probes reproduced findings 1–2.
- **9,500** deterministic Python checks against the written SET.1, OVL.2–5 and NEAR.2 formulas passed. No R used.
- No tracked cross-language harness is shipped; ignored development scratch files were not evaluated as HEAD artifacts.
- Worktree unchanged; no staged files. No independent R fidelity evaluation performed.
