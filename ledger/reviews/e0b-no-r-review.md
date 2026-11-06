# E0b-no-r review

HEAD: `7322367` (`exp/e0b-no-r`). Read-only review. No worktree edits or commits.

## Findings

Paths below are relative to `experiments/e0b-no-r/`.

1. **Major — Strand-aware intra-range operations lack `ignore_strand`.**
   `src/ranges_e0b/granges.py:284–296`
   TARGET requires strand-aware GRanges operations including `ignore.strand`. `resize`, `flank`, and `promoters` always apply the minus-strand transformation and expose no override. Probes of all three with `ignore_strand=True` raise `TypeError`. The API-presence tests do not check this required option.
   **Fix:** Add `ignore_strand=False`; bypass the minus-strand mask when true. Test both modes on minus-strand ranges against the existing IRanges kernels.

2. **Major — Mutable metadata breaks the documented immutable-value design.**
   `src/ranges_e0b/_sew.py:125–131`; `src/ranges_e0b/iranges.py:94–96`; `src/ranges_e0b/granges.py:172–174`
   Incoming Polars DataFrames are retained directly, and both `mcols` accessors return the internal object. Intra-range results share it. With `m = pl.DataFrame({'score': [1]})`, construct `x` using `mcols=m` and `y=x.shift(5)`: `m.replace_column(0, pl.Series('score', [99]))` changes both objects. Calling `x.mcols.insert_column(...)` also changes both. Metadata can change without creating a new range object, contrary to README:22–23.
   **Fix:** Clone caller-owned DataFrames on ingestion and return clones from public metadata accessors. Add constructor-alias and accessor-mutation regressions for both classes and derived results.

3. **Major — Replacing seqinfo cannot remove unused levels.**
   `src/ranges_e0b/granges.py:31–35` (called by `with_seqinfo`)
   `with_seqinfo` promises only that used seqnames must remain. Its validation accepts removal of unused levels, but `_remap` then looks up every old level, including unused ones. Reproduction: `GRanges('chr1', IRanges(start=1, width=3), seqinfo=Seqinfo(['chr1', 'chr2'])).with_seqinfo(Seqinfo(['chr1']))` raises `KeyError: "unknown seqlevel 'chr2'"`. Empty objects with old levels have the same problem.
   **Fix:** Remap only codes present in the object, with an early return for empty codes. Test unused-level removal, removal after subsetting, empty objects, and rejection of removed used levels.

4. **Minor — Empty Boolean masks bypass length validation.**
   `src/ranges_e0b/_sew.py:107–113`
   PLAN explicitly requires Boolean masks to match object length. The generic empty-array return precedes the Boolean shape check. For a two-range IRanges or GRanges, `x[np.array([], dtype=bool)]` silently returns zero ranges, while a one-element Boolean mask raises `ValueError`. An accidentally empty mask can silently discard data.
   **Fix:** Validate Boolean shape before handling generic empty indices. Preserve `x[[]]` as a valid empty positional subset. Test both classes.

## Completeness

- Every named TARGET operation has a public entry point; no whole operation is missing.
- Partial required support: intra-range `ignore_strand` (finding 1).
- Circular wraparound is explicitly unimplemented in README/NOTES. This is a recorded scope limitation, not a newly adjudicated R-fidelity finding.
- Public-API tests check names, not signatures or behavior completeness.

## Validation

From `experiments/e0b-no-r/`:

```sh
PYTHONDONTWRITEBYTECODE=1 HYPOTHESIS_STORAGE_DIRECTORY=/tmp/e0b-no-r-review-hypothesis UV_NO_SYNC=1 uv run pytest -p no:cacheprovider
```

**179 passed in 1.44 s.** Cache/bytecode writes were disabled or redirected outside the worktree. No cross-language harness ships with this experiment.

Foreground Python probes reproduced all four findings. An additional small exhaustive zero-width nearest probe matched the implementation's own documented distance rule; no finding raised.

Final `git status --porcelain` and `git diff --cached --name-only` were empty. No tests were added or updated.

## Approach and maintainability

- Visible artifacts fit the knowledge-only approach: independent numpy/Polars package, hand-derived examples, Python property-test oracles, and no R harness/goldens or ADR suite. No forbidden dependencies appear in project metadata. This does not audit commands from prior sessions.
- Scoped history shows staged commits with the required attribution. No other refs, worktrees, background material, R, or BiocPy were inspected.
- Shared array kernels and thin wrappers are a useful seam. Tests cover all operation groups and several invariants. Main gaps are metadata ownership, unused seqlevels, and argument-validation boundaries.
- No independent R-versus-Python fidelity evaluation was performed.

**Recommendation:** Fix findings 1–3 and add the boundary regressions before acceptance. Finding 4 is a small validation fix. No blocker findings.
