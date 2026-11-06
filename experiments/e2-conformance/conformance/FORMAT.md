# Conformance format

Normative for both codecs (`r/codec.R`, `src/ranges_e2/conformance/codec.py`).
JSON everywhere. R parses with `jsonlite::fromJSON(simplifyVector = FALSE)`,
Python with `json`. Both keep `1` (integer) distinct from `1.0` (double).
`null` = R `NA` (inside vectors) or `NULL` (elsewhere).
**All indices are 1-based** (R-native). Only `ops.py` converts for Python.

## Atomic vectors

Canonical (R-written) output is always tagged:

| Tag | R type | Notes |
|---|---|---|
| `{"$int": [1, null]}` | integer | |
| `{"$dbl": [1.5, "NaN"]}` | double | `"NaN"`, `"Inf"`, `"-Inf"` are strings |
| `{"$chr": ["a"]}` | character | |
| `{"$lgl": [true]}` | logical | |

Optional `"names": [...]` key on a tagged vector.

**Authoring shorthand (inputs only):** bare scalar/array. Integer literal →
integer; literal with `.` or exponent → double; string → character; bool →
logical; all-null → logical NA. Mixed kinds are an error. `[]` is ambiguous →
use a tag.

## Objects (`"$type"`)

| `$type` | Canonical fields | Notes |
|---|---|---|
| `IRanges` | `start` `width` (`$int`), `names` (`$chr`\|null), `mcols` (DFrame) | Inputs: any two of `start`/`end`/`width`. Any IntegerRanges subclass encodes as IRanges. |
| `GRanges` | `seqnames` (`$chr`), `start`, `width`, `strand` (`$chr` of `+ - *`), `names`, `mcols`, `seqinfo` (Seqinfo) | Level order lives in `seqinfo`. Inputs: `strand` default `*`; no `seqinfo` → levels in order of first appearance (decode rule, not R behavior). |
| `Seqinfo` | `seqnames` `$chr`, `seqlengths` `$int`, `isCircular` `$lgl`, `genome` `$chr` | NA → null |
| `DFrame` | `nrow`, `columns` {name: value} | Key order significant. `mcols(x) == NULL` encodes as `columns: {}`. |
| `Hits` | `from` `to` (`$int`), `nLnode`, `nRnode` | Any Hits subclass; pairs in R's emitted order |
| `Rle` | `values` (tagged or Factor), `lengths` `$int` | |
| `RleList` | `elements` {name: Rle} | |
| `IntegerList` | `elements` [[int,…],…], `names` | e.g. `revmap` |
| `Factor` | `levels` `$chr`, `codes` `$int` | |

Anything else an op returns is a harness error, not a pass.

## Cases

See `cases/*.json`. Per file: `group`, `clauses` {ID: statement},
`fixtures` {name: value}, `cases` [].

```json
{"id": "ir.reduce.mingap2", "clause": "IR.REDUCE.MINGAP", "op": "reduce",
 "args": [{"$ref": "x"}], "kwargs": {"min.gapwidth": 2}}
```

- `op`: R function name, run as `do.call(op, c(args, kwargs))` with
  GenomicRanges attached. Includes `"["`, `"c"`, constructors, replacement
  functions such as `"names<-"` (pass `value=` as kwarg).
- `kwargs` use R names. Python's adapter maps them; unknown kwargs error.
- `{"$ref": "name"}` → file-level fixture (fixtures may ref earlier ones).
- Optional `compare`: `exact` (default) or `hits-unordered` (set of pairs).
- Validator rules: unique ids across all groups; each case names exactly one
  declared clause; each clause has ≥1 case; refs resolve.

## Expected (`expected/<group>.json`, generated, never hand-edited)

```json
{"meta": {...}, "cases": {"<id>": {"value": <enc>, "warnings": ["…"]}}}
{"meta": {...}, "cases": {"<id>": {"error": {"message": "…", "class": ["…"]}}}}
```

- Error case passes in Python if the adapter raises `ValueError`, `TypeError`,
  `IndexError` or `KeyError`. Messages are not compared.
- Warnings are recorded, not compared. Python may warn freely.

## Comparison

Deep tree equality of encoded results (each side encodes its own result).
Numbers: `abs(a-b) <= 1e-9 * max(1, abs(b))`. Tags carry type; an empty
object equals an empty array.
