# E1 review — HEAD `43fcec8`

**Request changes.** No required operation family is wholly absent. Construction, metadata, named subsetting, flank, coverage, genomic overlaps/nearest, and parallel intersection remain partial. Existing artifacts identify 47 bug cases and five undecided divergences.

Paths below are relative to `experiments/e1-differential/`.

## Findings

1. **Major — `tests/differential/harness.py:108,117,139`: cache fingerprints never match.** Validation includes a script digest and differently spaced versions; persistence uses R's plain version string. Even a freshly refreshed cache fails validation. Every normal R-enabled session reruns all files. Offline mode also ignores script changes. **Fix:** persist and validate one fingerprint; check its script component offline. Add cache-hit and invalidation tests.

2. **Major — `src/ranges_e1/iranges.py:101,149`; `src/ranges_e1/_meta.py:117`: immutability is not enforced.** A writable base array changes stored coordinates through its read-only view. Caller-owned/accessor-returned Polars frames mutate both the original and shifted objects. Probe: one range acquires two metadata rows after `df.vstack(..., in_place=True)`. This contradicts ADR-0002. **Fix:** detach constructor buffers, protect accessor arrays, and clone metadata at public boundaries. Test mutation isolation.

3. **Major — `src/ranges_e1/granges.py:603`: circular overlap semantics are missing.** `_pairs` uses linear lifted coordinates without circularity. Existing `rt-circ-*` cases expose 17 failures across overlap types, count/any/subset, self-hits and nearest. **Fix:** handle periodic coordinates using merged circular seqinfo, preserving original indices and deduplicating hits; remove these xfails after passing.

4. **Major — `src/ranges_e1/_meta.py:120`: short metadata cannot recycle.** Twelve existing constructor/setter cases expect recycled rows, including nondivisor lengths; Python rejects them. **Fix:** implement the demonstrated row recycling and retain zero-row rejection on nonempty objects.

5. **Major — `src/ranges_e1/granges.py:94-99`: coordinate recycling leaves invalid names.** `rt-meta-gr-recycle-range-names` produces two ranges with one name. Independently verified: selecting its second row raises `IndexError`. **Fix:** recycle inherited range names alongside coordinates; assert names length and test subsequent subsetting.

6. **Major — `src/ranges_e1/iranges.py:507-514`; `src/ranges_e1/granges.py:799-803`: coverage input support is incomplete.** Existing cases show valid missing numeric weights rejected, IRanges weight-column lookup absent, and GRanges shift-column lookup absent. **Fix:** resolve metadata arguments before numeric normalization and preserve missing-weight effects locally in the run-length kernel, including overlapping ranges. Clear the seven related xfails.

7. **Major — `src/ranges_e1/granges.py:560`: existing `hit` metadata breaks `pintersect`.** Horizontal concatenation raises Polars `DuplicateError`; `rt-meta-gr-pintersect-existing-hit-keep` expects replacement in place. **Fix:** replace/add the column without duplicates while preserving other columns and the tested drop-nohit behavior.

8. **Major — `src/ranges_e1/iranges.py:264`; `src/ranges_e1/granges.py:349`: flank discards missing logical flags.** `np.asarray` strips masked-array masks. Both existing `rt-meta-*-flank-na-logical` cases silently calculate coordinates instead of rejecting NA. **Fix:** validate masks before conversion in a shared logical-vector helper.

9. **Major — `src/ranges_e1/granges.py:39`; `src/ranges_e1/seqinfo.py:59,67`: invalid inputs are silently reinterpreted.** Existing cases show literal strand `"NA"` becoming wildcard, numeric/text circularity becoming booleans, and numeric seqnames accepted. **Fix:** distinguish missing values from literal text and enforce character sequence names/logical circularity before coercion.

10. **Minor — `src/ranges_e1/iranges.py:29`; `src/ranges_e1/granges.py:57`: default empty construction is absent.** Both existing `rt-edge-*-default-empty` cases fail. **Fix:** support omitted arguments as empty objects without weakening incomplete nonempty-coordinate validation.

11. **Minor — `src/ranges_e1/_meta.py:198`: blank names incorrectly participate in lookup.** `rt-edge-ordering-subset-empty-name` selects rows where the committed reference errors. **Fix:** exclude empty names from valid lookup targets; preserve first-match behavior for valid duplicates.

12. **Major — `src/ranges_e1/granges.py:242,667-677`: five divergences lack decisions.** Existing seqinfo replacement and negative-coordinate nearest/precede/follow cases remain xfailed without ADR authorization. **Fix:** explicitly choose parity or document/pin accepted Python behavior through ADR-linked divergence cases. Do not blindly reproduce cross-sequence sentinel artifacts. Xfail alone does not close the R-normative resolution loop.

## Validation and process

- Cache-only suite: **3888 passed, 52 xfailed**, 10.04 s.
- Full fresh-R suite: **3888 passed, 52 xfailed**, 46.92 s. Redirected harness cache to `/tmp`; no worktree writes.
- Four existing cases rerun with `--runxfail`: **four failures**, confirming findings 3, 5, 7 and 10.
- Python-only probes confirmed finding 2 and the downstream subset crash in finding 5.
- Approach substantially followed: staged port, extensible YAML harness, live R reference, generated/adversarial cases, documented ADRs. Resolution remains unfinished. Four historical stage commits omit the requested attribution trailer.
- Design is compact: shared kernels, thin public wrappers, sparse coverage, reproducible generators. Cache lifecycle and shared mutable storage are the main additional maintenance defects.
- No independent R-fidelity hunt. R only reran shipped cases. No worktree edits, commits or staged files.
