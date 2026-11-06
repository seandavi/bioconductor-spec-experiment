# Result projection

One JSON shape per result kind, in **R conventions**: 1-based indices, NA as
`null`, typed vectors. R (oracle) and every adapter emit exactly this.
Shape is normative in `projection.schema.json`. Inputs handed to ports use the
same shapes.

## Scalars inside values

| R | JSON |
|---|---|
| integer | JSON integer |
| double | JSON number written `%.17g` (round-trip exact); NA → `null`; NaN → `"NaN"`; ±Inf → `"Inf"` / `"-Inf"` |
| character | string; NA → `null` |
| logical | `true` / `false`; NA → `null` |

`jsonlite::toJSON(digits = NA)` is **not** round-trip exact (15 digits);
the oracle must format doubles itself.

## Kinds

| `kind` | Fields | R source (oracle) |
|---|---|---|
| `vector` | `type` (`int` `dbl` `chr` `lgl`), `values`, `names` (`null` or list) | atomic vector; `names(x)` |
| `NULL` | — | `NULL` |
| `IRanges` | `start`, `width` (int lists), `names`, `mcols` | any `IntegerRanges` (IRanges, NormalIRanges, ...) |
| `GRanges` | `seqnames` (chr), `start`, `width`, `strand` (chr of `+ - *`), `names`, `mcols`, `seqinfo` | any `GenomicRanges` |
| `Seqinfo` | `seqnames` (chr), `seqlengths` (int\|null), `isCircular` (lgl\|null), `genome` (chr\|null) | `seqlevels`, `unname(seqlengths)`, `unname(isCircular)`, `unname(genome)` |
| `DataFrame` | `nrow`, `rownames` (`null` or chr), `columns` | `DataFrame` |
| `Factor` | `levels` (chr), `codes` (1-based int\|null) | `factor` (not inside an Rle) |
| `IntegerList` | `elements` (list of int lists), `names` | any `IntegerList` (e.g. `revmap`) |
| `Hits` | `from`, `to` (1-based int), `nLnode`, `nRnode` | any `Hits` / `SelfHits`; pairs in emitted order |
| `Rle` | `type` (`int` `dbl` `lgl`), `values`, `lengths` | non-factor `Rle` (coverage) |
| `RleList` | `names`, `elements` (list of `Rle` bodies: `type`, `values`, `lengths`) | `RleList` (GRanges coverage) |
| `error` | `message`; port side adds `type` (exception class) | R condition of class `error` |

## Rules

- `names`: `null` when the object has no names (R `NULL`). Otherwise a list the length of the object; NA names → `null` entries. A zero-length object always projects `names: null` (and `rownames: null`): `character(0)` vs `NULL` names carry no information.
- `mcols` / `columns`: ordered list of `[name, value]` pairs; value is a `vector` (with `names: null`), `Factor` or `IntegerList` projection. No columns → `[]` (R `NULL` mcols and 0-column DataFrame both).
- `end` is not projected: it is `start + width - 1` by definition.
- **Factor inside an Rle** (`seqnames(gr)`, `strand(gr)`) → `vector` of type `chr`, expanded per element. Rle is out of scope as a class (TARGET.md); levels live in `seqinfo`.
- `Rle` runs are **maximal**: no two adjacent runs with equal values (NA equals NA). Every length ≥ 1. A zero-length Rle has empty `values` and `lengths`.
- `type` is the R storage type. Adapters map the port's dtype faithfully (PROJECTION is not a place to coerce toward R): integer dtypes → `int`, floating → `dbl`, bool → `lgl`, str → `chr`.
- Index-valued results (`ops.json` `index_valued`) are 1-based. Adapters convert from the port's documented base.
- Hits from `findOverlaps(x)` (self) are `Hits` too; R's `SortedByQuery*` class is recorded by the oracle (`r_class`), not projected.
- Anything R returns outside this table → case `excluded: r_unprojectable:<class>`. Anything a port returns that the adapter cannot project → `adapter_error`.

## Oracle extras (expected.jsonl only)

| Field | Meaning |
|---|---|
| `r_class` | `class(result)[1]` |
| `sorted_by_query` | `is(result, "SortedByQueryHits")` (Hits only) |
| `universe` | Hits projection of the same call with `select = "all"` (only where COMPARATOR.md rule 3 applies) |
| `warnings` | R warning messages; recorded, never compared |
| `domain` | class of the first argument: `IRanges` \| `GRanges` \| `Seqinfo` \| `other` |
| `features` | derived from canonical inputs: `empty`, `zero_width`, `names`, `mcols`, `strand`, `multi_seqlevel`, `seqlengths`, `circular`, `na_arg`, `r_error` |
