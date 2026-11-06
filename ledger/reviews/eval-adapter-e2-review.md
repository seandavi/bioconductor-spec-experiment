# E2 adapter review

## Thinness

1. **High — `eval/adapters/e2.py:224`**: Every list column named `revmap` is incremented, including ordinary input metadata and results of non-index-valued operations. Probe: `project(build_IRanges(p))` changes user `revmap=[[5]]` to `[[6]]`. **Fix:** distinguish port-generated index revmaps from ordinary IntegerList metadata; shift only the former, preserving provenance through subsequent calls.

2. **Medium — `eval/adapters/e2.py:276,432,514`**: `dict(zip(names, values))` silently deduplicates named constructor seqlengths, gaps bounds and coverage widths. Probes with duplicate `chr1` names retain only the last value. **Fix:** declare duplicate-name mapping forms Unsupported before conversion; do not silently discard entries.

3. **Medium — `eval/adapters/e2.py:409-410`**: `flank` accepts but completely ignores `use_names`. Probe with `FALSE` still returns the input names. Neither public flank method accepts this option. **Fix:** reject supplied `use_names` with `UnknownArgument`, consistent with the other unsupported `use.names` forms; do not implement name removal in the adapter.

No exception substitution, output sorting, clipping or range re-derivation found.

## Translation

4. **Medium — `eval/adapters/e2.py:73-75,516-522`**: Numeric-list conversion loses explicit R storage type. For valid `RVector('dbl', [1])`, coverage emits integer Rle; the same port call with a float64 array emits double Rle. **Fix:** preserve declared numeric dtype when translating coverage weights, including integral-valued doubles and empty vectors; retain NA rather than coercing it away.

5. **Medium — `eval/adapters/e2.py:340-347`**: `seqlengths`, `isCircular` and `genome` drop the sequence names from their vector projections. These are named R accessors; the port's Seqinfo table supplies the matching names. **Fix:** supply `names=list(si.seqnames)` to each returned RVector. Do not name `seqlevels` itself.

6. **Medium — `eval/adapters/e2.py:302-303`**: Supported DataFrame inputs/mcols results use `_DF`, which has no `__len__`; `length(DataFrame)` raises an adapter TypeError. **Fix:** translate this case to the public polars frame's column count (`frame.width`), not its row count.

## Unsupported declarations

7. **Medium — `eval/adapters/e2.py:305-307,389-391`**: Seqinfo `seqnames` and character subsetting are unnecessarily Unsupported. Public `Seqinfo.seqnames` and `Seqinfo.subset(names)` provide both, including unknown-name insertion. Confirmed with direct probes; resolved inputs include `si.acc.seqnames`, `si.subset.names` and `si.subset.unknown-name-adds-level`. **Fix:** add Seqinfo dispatch to those public APIs, projecting seqnames as a character vector.

8. **Medium — `eval/adapters/e2.py:112-115`**: Constructor metadata rejects native Factor arguments, although the port accepts polars Enum columns and the adapter already builds/projects that representation elsewhere. Probe: Factor metadata is Unsupported through `op_IRanges`, but the equivalent public constructor with an Enum column succeeds. **Fix:** translate Factor metadata to Enum using the existing level/code conversion, shared with input construction.

9. **Low — `eval/adapters/e2.py:273-275`**: Explicit `seqlengths=NULL` is declared Unsupported, although public `GRanges(..., seqlengths=None)` accepts it. **Fix:** pass None through; retain the unsupported declaration for unnamed non-NULL vectors.

## Validation

- `cd eval && just roundtrip e2`: exit 0; **5866 distinct inputs**, no reported failures. The additional user-revmap probe above fails roundtrip outside that corpus.
- `cd eval && just adapter-check e2`: exit 0; **9650 cases: 7359 value, 1570 port_error, 721 unsupported; 0 adapter_error**. No scored results or expected values read.
- `uv run --frozen --project /Users/davsean/Documents/git/bioconductor-spec-experiment.worktrees/e2-conformance/experiments/e2-conformance python /tmp/e2-review-probe.py`: exit 0; targeted adapter/public-API probes confirmed the examples above. Probe log: `/tmp/e2-review-probe.log`.
- Final `git status --porcelain=v1` and `git diff --cached --name-only`: empty. No source edits, commits or staged files. Required dev-loop commands generated only ignored runner artifacts.
