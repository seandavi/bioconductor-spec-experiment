# Case schema

Shape is normative in `case.schema.json`; this file gives the semantics.
A case holds **inputs and a call only**. No expected values, ever.

## Files

| File | Written by | Read by |
|---|---|---|
| `eval/cases/{e1,e2,e4,heldout}.jsonl` | converters / heldout generator (committed) | oracle |
| `eval/cases/unconvertible.jsonl` | converters (committed) | report |
| `build/resolved.jsonl` | oracle | port runner (the only case file a port sees) |
| `build/expected.jsonl` | oracle | comparator only |

One JSON object per line. Lines sorted by `id`.

## Case

| Field | Meaning |
|---|---|
| `id` | `<suite>/<file-stem>/<origin id>`; heldout: `heldout/<op>/<k>`. Unique. |
| `source` | `e1` \| `e1-redteam` \| `e2` \| `e4` \| `heldout` |
| `origin` | `{suite, file, id, tags}`: where it came from (file relative to the suite root) |
| `op` | top-level `call.op` (indexes `ops.json` → family) |
| `inputs` | ordered `{name: Literal}`; a literal may `$ref` an earlier input |
| `call` | `Call` (below) |

## Literals (inputs only)

R builds every input with the fixed rule below (the "eval decode"), then
projects it (PROJECTION.md). Ports only ever see the projection.

| Literal | R decode |
|---|---|
| `{"int"\|"dbl"\|"chr"\|"lgl": [..], "names"?: [..]}` | typed atomic vector; `null` element = NA; dbl also takes `"NaN"`, `"Inf"`, `"-Inf"` |
| `{"IRanges": {start?, end?, width?, names?, mcols?}}` | `IRanges(start=, end=, width=, names=)` (absent = not passed), then `mcols(x) <- DataFrame(cols, check.names=FALSE)` if `mcols` |
| `{"GRanges": {seqnames, start?, end?, width?, names?, strand?, seqinfo?, mcols?}}` | `GRanges(seqnames=, ranges=IRanges(start,end,width,names), strand= (if given), seqinfo= (if given))`, then mcols as above |
| `{"Seqinfo": {seqnames, seqlengths?, isCircular?, genome?}}` | `Seqinfo(...)` with the given args |
| `{"DataFrame": {columns, nrow?}}` | no `nrow`: `DataFrame(cols, check.names=FALSE)`; with `nrow`: `make_zero_col_DFrame(nrow)`, then columns added in order |
| `{"Factor": {levels, codes}}` | `factor(levels[codes], levels=levels)`; codes 1-based, `null` = NA |
| `{"IntegerList": {elements, names?}}` | `IntegerList(list(...))` |
| `{"$ref": name}` | an earlier input (e.g. a shared Seqinfo) |

- `mcols` / `columns`: ordered list of `[name, value]` pairs; value is a vector, Factor or IntegerList literal.
- Fields inside object literals are vector literals (a bare JSON array is not allowed: type must be explicit).
- An input R cannot build makes the case `excluded: r_input_error` (counted, listed). Construction itself is tested via constructor **ops**, never via inputs.

## Call

```json
{"op": "findOverlaps", "args": [{"$ref": "q"}, {"$ref": "s"}], "kwargs": {"type": "within", "maxgap": 2}}
```

`args` / `kwargs` values (`Arg`):

| Arg | Meaning in R | Port runner passes |
|---|---|---|
| `{"$ref": name}` | the input | the adapter-built native object |
| `{"$call": Call}` | nested call result | the adapter's native result |
| vector literal | typed vector | `RVector` |
| `"s"` / `true` / `3` / `2.5` | length-1 chr / lgl / int / dbl | `RVector` of length 1 |
| `null` | `NULL` | `None` |

- Bare numbers: integral → int (must fit int32); non-integral → dbl. Integral doubles and out-of-int32 values are written `{"dbl": [..]}`. (JSON Schema cannot tell `2` from `2.0`; `spec/examples/validate.py` checks it.)
- `kwargs` keys are R argument names (`ignore.strand`, `min.gapwidth`); order kept.
- R evaluates `do.call(op, c(args, kwargs))` with IRanges + GenomicRanges attached.
- Object literals never appear in calls: converters hoist them into `inputs` (`_a1`, `_a2`, ...).

## Hash

`sha256` of canonical JSON (`sort_keys`, no whitespace) of `{"inputs", "call"}`.
Ids, tags and source are not hashed. The oracle cache keys on
(R + package versions, oracle script digest, case hash).
