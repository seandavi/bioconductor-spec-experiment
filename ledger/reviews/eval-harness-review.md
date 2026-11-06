# Oracle review

**Gate: changes requested.** Three reproducible defects; no worktree edits or commits. Reviewed `ef16e77`. Findings below use synthetic probes, not evidence that existing committed scores are already affected.

## Findings

1. **High — projection failure can earn a both-error pass.** `eval/src/bioc_eval/compare.py:163–164`; classification at `eval/runner/run_port.py:69–74`.
   - Live R rejects `IRanges(start=1, width=-1)`. A probe adapter whose operation **returns successfully**, then calls a failing native getter during projection, produces `port_error`, `phase=project`; comparator returns **pass/R1**. A projection timeout takes the same path. Thus broken translation/projection can conceal a port accepting an invalid call.
   - **Fix:** restrict R1 to operation-evaluation errors (`phase=call`). Projection failures must never earn R1; classify adapter projection failures as `adapter_error`, retaining a distinct unexpected-error classification for genuine native getter failures if desired. Update COMPARATOR/ADAPTERS and add projection-exception/timeout regressions against an R-error case.

2. **Medium — case hashes erase semantically significant argument order.** `eval/src/bioc_eval/common.py:22–29`; deduplication at `eval/src/bioc_eval/oracle.py:197–211`.
   - `IRanges(start=1,width=2,aa=3,bb=4)` and the same kwargs ordered `bb,aa` have identical `case_hash`. Live R returns metadata columns `[aa,bb]` versus `[bb,aa]`. Deduplication/cache reuse therefore supplies one call's expected projection to the other, although column order is explicitly exact.
   - **Fix:** preserve ordered `inputs`/`kwargs` recursively in the hash representation (e.g. ordered pairs), while canonicalizing only unordered structural mappings. Update CASES.md and invalidate old keys. Regression: execute both cases in one oracle build and again from cache; each must retain its own column order.

3. **Medium — R JSON projection deletes valid string content.** `eval/oracle/oracle.R:13–20`, especially line 19.
   - Live runner changes the range name `"a\bb\fc\u0001d"` into `"abcd"`. Control characters other than newline/CR/tab are removed rather than escaped. This corrupts both resolved inputs and expected values, masking string-handling defects and making correct character results mismatch literal call arguments.
   - **Fix:** JSON-escape every control character, without deletion; using jsonlite for scalar strings does not compromise the custom exact-double writer. Add names and character-mcols round-trip tests covering backspace, form feed and U+0001–U+001F.

## Checks and evidence

| Area | Result |
|---|---|
| Converter fidelity | Replayed **all 6,140 converted cases** through their original R input/call decoders and eval's live R decoder/projection: **0 result/error-occurrence differences**. No experiment expected values read; warnings/error messages not compared. Eleven unconvertibles are listed separately. |
| Heldout | Live R ran **3,510 cases**, **0 exclusions**, 440 normative errors. Every TARGET op family represented; empty, zero-width, names, mcols, strands and seqinfo covered. |
| Projection completeness | Recursive live outputs/inputs exercised all 12 declared kinds, including nested Factor/IntegerList; exception: string fidelity defect above. |
| Comparator | Independent value/name/order/NA/storage-type mutations all mismatch. Focused tests check R2 preserves query order, multiplicity and dimensions; R3 rejects out-of-universe picks and wrong names/types/lengths. R2's within-query relaxation is narrow; live R also confirms the documented nearest/all inconsistency motivating R3's own-pick exception. |
| Self-tests | On freshly generated live-R heldout output: **r_echo 3,510/3,510 pass**; mutant **0 prediction discrepancies**, including both relaxation controls. These tests do not detect finding 1. |
| Adapter instructions | Sufficient for a fresh author: signatures, RVector/NULL/missing distinctions, indexing, projection helpers, unsupported forms, forbidden semantic repairs, dev loop. Amend classification for finding 1 before delegation. No real port adapters reviewed. |

Commands: `PYTHONDONTWRITEBYTECODE=1 eval/.venv/bin/python -m pytest -p no:cacheprovider eval/tests/{test_compare,test_binding,test_classification,test_helpers,test_converters}.py -q` (**69 passed**); `Rscript /tmp/eval-oracle-review/converter-probe.R`; direct `Rscript eval/oracle/oracle.R` heldout/edge runs; `run_port.py --port r_echo` and `--port mutant` with all outputs under `/tmp/eval-oracle-review`; `edge-probes.py` there reproduces all three findings. Probe scripts/results remain there for independent reproduction.

R: 4.6.0; IRanges 2.47.2; GenomicRanges 1.65.1. `git status --short` and `git diff --cached --name-only` were empty. No frozen port environments changed. Generator-author blindness and post-generation tuning caveats are disclosed in NOTES.md; historical blindness was not independently verified.

**Next:** fix the three findings, rerun regressions and self-tests, then obtain reviewer approval before scoring real ports.

Co-Authored-By: Claude <noreply@anthropic.com>
