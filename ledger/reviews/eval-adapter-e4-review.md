## Thinness
No forbidden semantic shims found. Index-base shifts, Factor encoding and Rle normalization are permitted representation changes; no adapter sorting, deduplication, clipping or exception-to-value substitution found.

## Findings
1. **Medium — translation:** `eval/adapters/e4.py:452` discards named coverage widths. Existing input-only case `e2/granges-coverage/gr.cov.width.named-reordered` produces lengths `{chr1:5, chr2:3, chrX:12, chrU:6}` through the adapter, versus `{chr1:6, chr2:12, chrX:3, chrU:5}` using the port's documented `GRanges.coverage(width=dict)` API. **Fix:** translate named width RVectors to the public seqname→width mapping; leave unnamed widths positional.
2. **Medium — translation/projection:** `eval/adapters/e4.py:212-213` ignores DataFrame `nrow` and `rownames`. Probe: a valid three-row, zero-column input projects back with `nrow:0`. Non-null row names are also silently dropped. **Fix:** preserve requested height using public Polars construction (a temporary column followed by `drop` preserves zero-column height); explicitly reject non-null rownames if the native representation cannot retain them. Existing resolved inputs do not exercise these forms.
3. **Medium — false unsupported:** `eval/adapters/e4.py:253-255` rejects `length` on exported port `Rle` and `Hits`, although both provide public `__len__`. Existing case `e1/coverage/cov-length` is rejected despite its native Rle having length 14; a Hits probe similarly rejects native length 1. **Fix:** admit these types and delegate to `len(x)`. The five static seqinfo-setter unsupported declarations agree with the inspected public API.

## Validation
- `cd eval && just roundtrip e4`: passed, 5,866 distinct inputs.
- `cd eval && just adapter-check e4`: passed; 9,650 cases, 7,214 values, 1,273 port errors, 1,163 unsupported, zero adapter errors. This is not a scoring run.
- `/tmp/e4_adapter_review_probe.py`: confirmed all three findings against public APIs.
- `git status --short`, `git diff --exit-code`, `git diff --cached --exit-code`: clean; no tracked edits, staged files or commits. No expected/scored results, port tests or goldens read.
