# Differential tests: how to add a case

A case runs in **R** (normative) and **Python**; canonical results are compared.
R results are cached in `r_cache/<file>.json` (committed). Never hand-edit the cache.

## Workflow
1. Append a case to `cases/<group>.yaml` (or a new `cases/*.yaml`).
2. `uv run pytest tests/differential -k <id>` — runs R for new/changed cases, writes the cache.
3. Fail? `uv run python tests/differential/show.py <id>` prints R vs Python + first difference.
4. Fix Python, or write `adrs/NNNN-*.md` and add `divergence:` (below). Commit case + cache.

`--refresh-r` reruns everything in R. R version change → all entries rerun → `git diff r_cache/`.
`RANGES_E1_NO_R=1` forces cache-only mode (stale cases skip).

## Case syntax
```yaml
cases:
  - id: reduce-zero-width        # unique across all files
    tags: [reduce, zero-width]
    inputs:
      x: {IRanges: {start: [5, 1, 20], width: [0, 10, 0]}}
    call: {fn: reduce, args: [$x], kwargs: {min.gapwidth: 2}}
```
- `fn`, kwargs use **R names**. `$x` = input; `{call: {...}}` nests; `$$` escapes a literal `$`.
- Literals: `{IRanges: {start,end,width,names,mcols}}`,
  `{GRanges: {seqnames,start,end,width,strand,names,mcols,seqinfo}}`,
  `{Seqinfo: {seqnames,seqlengths,is_circular,genome}}`, `{DataFrame: {col: [..]}}`,
  scalars, arrays, `null` = NA. Empty/ambiguous vectors: `{int: []}`, `{num: [..]}`, `{chr: []}`, `{lgl: []}`.
  Whole-number JSON numbers are R integers; use `{num: [..]}` for doubles.
- **Indices are R's (1-based) in cases and in canonical results.** The adapter converts.
- Escape hatch: `raw: {r: "<R expr>", py: "<Python expr>"}`; inputs are bound as variables
  (`np`, `pl` and `ranges_e1` names too). Python-side indices are 0-based there.
- Optional keys: `check_warnings: true` (both sides must warn, or both not);
  `compare: arbitrary` + `universe: {fn: findOverlaps, ... select: all}` (Python pick must be
  a member of R's hit set; needs an ADR); `divergence: {adr: "0004", python: <canonical>}`
  (asserts Python equals the ADR-sanctioned value; `{error: true}` for "Python errors").

## Comparison
Structural; int vs num by value (1e-9), everything else strict. Errors: pass iff both error
(messages ignored). Warnings ignored unless `check_warnings`. mcols column order matters.

## Canonical values (what both sides emit)
`IRanges {start,width,names,mcols}`; `GRanges {seqnames,start,width,strand,names,mcols,seqinfo}`;
`Hits {from,to,nLnode,nRnode}`; `Rle {values,lengths}`; `RleList`; `List {elements}`;
`vector {type: int|num|lgl|chr, values, names}`. A factor/factor-Rle (seqnames, strand) is a chr vector.
mcols column `revmap` is index-valued (+1 on the Python side).

## New R function in the adapter
Usually nothing: `fooBar(x, a.b=1)` → `x.foo_bar(a_b=1)`; `fn<-` → `set_fn` (value via `kwargs: {value: ..}`);
properties are read as attributes. Add an `ADAPTERS` entry in `adapter.py` only for
index-valued results (`post=_post_idx`), non-method functions, or argument reshaping.
Unit-test hooks in `tests/unit/test_adapter.py`.

## Files
`harness.py` (cases, cache, R runner, compare) · `adapter.py` (call table, encoder) ·
`codec.py` (literal decode, diff) · `r/run_cases.R` + `r/encode.R` (R side) · `show.py`.
