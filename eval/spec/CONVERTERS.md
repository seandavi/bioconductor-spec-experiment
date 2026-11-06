# Converters (E1, E2, E4 → eval cases)

Converted cases keep **inputs and calls only**. Each converter reproduces
its source's own R decoding as explicit eval literals, so R evaluates the
same construction the source suite did. Anything that cannot be expressed
becomes one line in `eval/cases/unconvertible.jsonl`:
`{source, origin, reason, detail}`. Nothing is dropped silently.

Suites are read from the frozen tags (READ-ONLY):

| Suite | Root (under `<worktrees>/<dir>/experiments/<dir>/`) | Files | Cases at freeze |
|---|---|---|---|
| E1 | `tests/differential/cases/` | `*.yaml` (PyYAML `safe_load`, as E1's harness) | 4515 (119 tagged `redteam`, 4 `raw`) |
| E2 | `conformance/cases/` | `*.json` | 921 |
| E4 | `cases/{iranges,granges}/` | `*.yaml` (PyYAML `safe_load`; E4's R reads `y`/`n` as strings too) | 715 |

## Common rules

| Rule | Detail |
|---|---|
| id | `<suite>/<file stem>/<origin id>` (E4 stem includes domain: `e4/granges-overlaps/<id>`) |
| dropped fields | everything expected-adjacent: E1 `divergence`, `check_warnings`; E2 `compare`, `xfail.json`, `expected/`; E4 `python`, `adr`, `compare`, `oneof_from`; all goldens/caches |
| tags | `origin.tags`: E1 `tags`; E2 `[group, clause]`; E4 `[domain, clause]` |
| vocabulary | every op in the call tree must be in `ops.json`, else `op-out-of-vocab:<op>` (e.g. `runValue`, `runLength`, `as.vector`, `[[`) |
| inputs | only inputs the call references (transitively), in source order |
| hoisting | object literals inline in `args`/`kwargs` → inputs `_a1`, `_a2`, ... → `$ref` |
| identity rewrite | `identity($x)` with `x` an IRanges/GRanges/Seqinfo literal (E1 `grid_construct`, E4 construct cases) → constructor-op tree mirroring the source's decoder (below), so construction errors are still tested. `mcols` → `mcols<-`(..., value = `$ref` hoisted DataFrame). `identity` on anything else → `op-out-of-vocab:identity` |
| scalars | int in int32 → bare int; int outside int32 → `{"dbl":[v]}`; integral float → `{"dbl":[v]}`; other float → bare; str, bool → bare |

Constructor trees (identity rewrite):

```
IRanges literal → IRanges(start=, end=, width=, names=)            # only given fields
GRanges literal → GRanges(seqnames=, ranges=IRanges(...), strand=, seqinfo=Seqinfo(...))  # strand/seqinfo only if the source decoder passes them
+ mcols         → `mcols<-`(<above>, value = $ref _dfN)
```

## E1

| Source form | Eval |
|---|---|
| `source` | `e1-redteam` if file stem starts `redteam_` or tags contain `redteam`; else `e1` |
| `call: {fn, args, kwargs}` | `{op: fn, args, kwargs}` |
| `"$x"` / `"$$x"` / other string | `{"$ref": "x"}` / chr `"$x"` / chr scalar |
| `{call: {...}}` | `{"$call": {...}}` |
| `null` | `{"lgl": [null]}` (E1 `decode(NULL)` is `NA`) |
| `{NULL: ...}` | `null` |
| `{int\|num\|chr\|lgl: [..]}` | `{int\|dbl\|chr\|lgl: [..]}` |
| untyped list | R `unlist`: all bool → lgl; all int32 ints → int; ints + floats or out-of-int32 → dbl; all str → chr; all null → lgl; `[]` → int (E1 default); other mixes → `literal:mixed` |
| `IRanges{start,end,width,names,mcols}` | same fields, each decoded as above (no coercion: E1 passes types through) |
| `GRanges{...,strand?,seqinfo?}` | same; `strand`/`seqinfo` only if present |
| `Seqinfo{seqnames, seqlengths, is_circular, genome}` | `isCircular` renamed; `seqlengths` → int (E1 applies `as.integer`; non-numeric → `literal:coerce`) |
| `DataFrame{col: v}` / `mcols{col: v}` | `columns` pairs, no `nrow` |
| `raw` (4) | hand-mapped table below; any other `raw` → `raw-r` |

| E1 raw id | R text | Eval call |
|---|---|---|
| `rt-edge-ir-default-empty` | `IRanges()` | `IRanges()` |
| `rt-edge-gr-default-empty` | `GRanges()` | `GRanges()` |
| `rt-meta-seqinfo-numeric-names` | `Seqinfo(c(1L, 2L))` | `Seqinfo({"int":[1,2]})` |
| `rt-meta-gr-recycle-range-names` | `GRanges(c("a","b"), IRanges(1L, width=2L, names="only"))` | `GRanges({"chr":["a","b"]}, $call IRanges(1, width=2, names="only"))` |

## E2

| Source form | Eval |
|---|---|
| case `{id, clause, op, args, kwargs}` | `{op, args, kwargs}`; `source: e2` |
| `{"$ref": f}` | fixture `f` copied into `inputs` (with its transitive refs, fixture order) |
| `{"$int"\|"$dbl"\|"$chr"\|"$lgl": [..], names?}` | typed vector (`$dbl` keeps `"NaN"`/`"Inf"`) |
| shorthand array/scalar | E2 rules: all int → int; any float → dbl; str → chr; bool → lgl; all null → lgl NA; `[]` → `literal:ambiguous-empty` (E2 itself rejects it); mixes → `literal:mixed` |
| `$type: IRanges` | IRanges literal; `mcols` DFrame → `columns` |
| `$type: GRanges` | GRanges literal with E2's decode made explicit: no `seqinfo` → `Seqinfo(unique(seqnames))` in appearance order; no `strand` → `rep("*", n)` |
| `$type: Seqinfo` / `DFrame` / `Factor` / `IntegerList` | Seqinfo / DataFrame (`nrow` kept) / Factor / IntegerList literal |
| `$type: Hits` / `Rle` / `RleList` as input | `input-kind:<k>` |

## E4

| Source form | Eval |
|---|---|
| `call: {op, args: [names], kwargs}` | args → `$ref`; `source: e4` |
| kwarg `{ref: name}` / list / scalar / `null` | `$ref` / untyped vector (E4 `vec`: as E1 but `[]` → lgl) / scalar / `null` (R `NULL`) |
| `IRanges{start,end,width}` | **int** (E4 applies `as.integer`: floats truncate, bool → 0/1, null → NA, str → `literal:coerce`); `names` → chr; `mcols` untyped |
| `GRanges{seqnames, strand?, seqinfo?, ...}` | seqnames/strand chr; `seqinfo: <string>` → `{"$ref": name}`; `strand`/`seqinfo` only if present |
| `Seqinfo{seqnames, seqlengths, isCircular, genome}` | chr / int / lgl / chr |
| `int\|num\|chr\|lgl` input | typed vector |

## Outputs

`eval/cases/{e1,e2,e4}.jsonl`, `eval/cases/unconvertible.jsonl`,
`eval/cases/conversion-summary.json` (per source: read, converted,
unconvertible by reason). All schema-validated (`case.schema.json`).
