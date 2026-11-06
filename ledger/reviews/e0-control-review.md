## Review

Reviewed `exp/e0-control` at HEAD `524cb18`. The supervisor confirmed a clean worktree and ran `uv run pytest -q`: **300 passed in 2.19s**. No files were edited. No independent R-fidelity investigation was performed.

### Correct

- **Completeness:** No required operation is missing from the public API. Construction, metadata, subsetting, intra/inter-range operations, set operations, overlaps, nearest, coverage, ordering, and Seqinfo are represented. Some behavior is partial or incorrect as detailed below.
- **Approach:** The inspected artifacts follow E0: ordinary pytest examples and Python-only property tests, with no cross-language harness, shared cases, or R-generated goldens. No prohibited imports were found. The supplied commit history shows staged implementation.
- **Maintainability:** Private array kernels and the GRanges grouping seam provide useful separation. README documents conventions, deliberate deviations, and maintenance. Validation is less consistently centralized than the algorithms.

### Numbered findings

1. **Major (P1): Empty ranges can select themselves as their nearest neighbor.**  
   **Locations:** `experiments/e0-control/src/ranges_e0/_nearest.py:84–94,108–109`; `experiments/e0-control/src/ranges_e0/granges.py:780–789`.

   Self-exclusion only filters the initial overlap results. Directional fallback can then select the same zero-width range. Supervisor-confirmed reproduction:
   ```python
   x = IRanges(5, width=0)
   x.nearest().tolist()                         # [0]
   x.nearest(select="all").to_numpy().tolist()  # [[0, 0], [0, 0]]
   ```
   The GRanges equivalent also returns itself. This contradicts the behavior asserted by `tests/test_iranges_nearest.py:97–101` and additionally produces duplicate IRanges hits.

   **Fix:** Exclude original query identity during directional candidate selection, before selecting the nearest candidate. Deduplicate combined pairs. Add singleton-empty and mixed empty/non-empty self-query regressions for both classes.

2. **Major (P1): GRanges distance silently discards rows for incompatible lengths.**  
   **Location:** `experiments/e0-control/src/ranges_e0/granges.py:816–822`.

   Any shorter operand is replaced with repeated copies of its first row, even when its length exceeds one. The subsequent length check therefore cannot reject incompatible inputs. This contradicts README line 100 and the method’s stated length-one recycling convention.

   Confirmed reproduction:
   ```python
   a = GRanges("chr1", IRanges([1, 100], width=1))
   b = GRanges("chr1", IRanges([10, 20, 30], width=1))
   a.distance(b).tolist()  # [8, 18, 28], silently ignores a's second row
   ```

   **Fix:** Validate original nonempty lengths before broadcasting: each must equal one or the common output length. Add unequal-length tests in both argument orders.

3. **Minor (P2): Coverage weights are not consistently recycled or validated.**  
   **Locations:** `experiments/e0-control/src/ranges_e0/granges.py:837–849`; `experiments/e0-control/src/ranges_e0/_coverage.py:35–39`.

   For two ranges, `g.coverage(weight=[2])` raises `IndexError`, despite the documented length-one recycling convention. Conversely, `g.coverage(weight=[1, 2, 999])` silently ignores the extra weight: GRanges subsets weights before the kernel can validate their original length. The IRanges kernel also rejects a length-one weight vector.

   **Fix:** Normalize scalar/length-one weights and validate their dimensionality and length before grouping. Share this normalization between IRanges and GRanges. Test scalar, length-one, per-range, undersized, and oversized weights.

4. **Minor (P2): Flank accepts incompatible argument lengths and creates malformed objects.**  
   **Locations:** `experiments/e0-control/src/ranges_e0/iranges.py:154–157`; `experiments/e0-control/src/ranges_e0/granges.py:346–349`.

   Confirmed reproduction:
   ```python
   y = IRanges(10, 20).flank(2, start=[True, False])
   len(y), y.start.tolist(), y.width.tolist()
   # (2, [8, 21], [2])
   ```
   NumPy broadcasting expands starts without expanding widths or associated metadata. The GRanges equivalent has the same inconsistent range arrays. This violates the vector-length invariant rather than rejecting an unsupported argument length.

   **Fix:** Validate `start` as one-dimensional and recycle only scalar/length-one inputs; otherwise require exactly the input range count. Use a shared boolean validator and add rejection tests for both classes.

### Merge verdict

**OK with notes — public scope and existing tests are complete, but the major correctness findings should be fixed before release.**