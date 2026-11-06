## Thinness

No score-inflating repairs found: no sorting, deduplication, clipping, recomputation or exception substitution. Coverage encoding, index shifts, combine→range and SeqInfo getter/setter composition are representation/public-API translations.

## Findings

1. **Medium — `eval/adapters/biocpy.py:282,297`**: Constructor metadata loses atomic dtype by converting every column to a Python list. An empty integer column projects as `IntegerList`; an all-NA integer column raises projection `ValueError`. **Fix:** translate `RVector` columns into typed/masked arrays using the existing `_column` conversion; pass native port columns through unchanged. Do not infer the original dtype from values.

2. **Medium — `eval/adapters/biocpy.py:94,129,275,277`**: Nested native results are rejected by translations that assume `RVector`. `IRanges(start=start(x), width=width(x))` raises `AttributeError: ndarray has no attribute list`; `x[order(x)]` reaches the non-RVector `Unsupported("subscript:NULL")` branch. **Fix:** accept native typed/masked vectors as well as literal `RVector` arguments, retaining NA and applying exactly one subscript base shift.

3. **Medium — `eval/adapters/biocpy.py:201`**: `IntegerList` input elements containing NA are passed to `np.array(..., dtype=int32)`, raising adapter `TypeError` before the port receives them. **Fix:** construct masked integer arrays for missing elements, then retain those masks through `p_integerlist` projection.

4. **Medium — `eval/adapters/biocpy.py:186`**: Blanket `seqlevels<-` unsupported declaration hides the public `SeqInfo.set_seqnames` API. At least unnamed positional renaming is directly available (`['a','b']`→`['c','d']`); GRanges can compose it with `get_seqinfo`/`set_seqinfo`, like the existing column setters. **Fix:** expose positional renaming; restrict unsupported declarations to argument forms the public API cannot represent. Do not implement name matching, dropping or adding levels yourself.

5. **Medium — `eval/adapters/biocpy.py:190-191`**: Public dependency APIs `biocutils.duplicated` and `biocutils.unique` accept both IRanges and GenomicRanges, despite the blanket unsupported declarations. Probes execute successfully (and retain duplicates incorrectly). **Fix:** dispatch directly to these public generics and expose their actual results; do not repair their equality/deduplication behavior.

6. **Low — `eval/adapters/biocpy.py:289,301`**: Empty `GRanges()` and `Seqinfo()` are declared unsupported although documented public `GenomicRanges.empty()` and `SeqInfo.empty()` provide them. **Fix:** route the no-argument forms to these factories, without supplying R defaults to other constructor forms.

## Validation

- `cd eval && just roundtrip biocpy`: **failed**, exit 1; 5,866 distinct inputs, 10 failures: four NA names become string `"None"`, six int32-end construction rejections. These are observed port limitations; do not repair them in the adapter.
- `cd eval && just adapter-check biocpy`: **passed**, exit 0; 9,650 cases: 6,354 value, 1,690 unsupported, 1,606 port_error, **zero adapter_error**. Findings 1–3 were exposed by additional probes, not this corpus.
- `PYTHONPATH=eval/runner:eval/adapters eval/ports/biocpy/.venv/bin/python /tmp/biocpy-review-probes.py`: completed; reproduced constructor dtype/NA failures, native constructor failure, IntegerList NA failure, public factories/rename and dependency generics. Logs: `/tmp/biocpy-review-{roundtrip,adapter-check,probes}.log`.
- Worktree source unchanged; no commits or staged files. Adapter-check generated its normal ignored `eval/build/results/biocpy.jsonl`. No `just run`, `just report`, expected values, port tests or other adapters were read/run.
