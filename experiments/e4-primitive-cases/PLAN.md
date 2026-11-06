# PLAN — e4 primitive shared cases

Port IRanges/GRanges (TARGET.md) to `ranges_e4`. Conformance lives in
YAML case files of primitives; R computes goldens; Python checks them.
Objects are never serialized — only constructor arguments in, primitive
projections out.

## 1. Layout

```
experiments/e4-primitive-cases/
  PLAN.md  README.md  NOTES.md  clauses.md  justfile
  pyproject.toml  uv.lock
  adr/NNNN-slug.md            # short decisions (non-exact compare, divergence)
  cases/iranges/<group>.yaml  # hand-authored: inputs, call, compare
  cases/granges/<group>.yaml
  goldens/iranges/<group>.json  # R-generated; never hand-edited
  goldens/granges/<group>.json
  runners/r/run_cases.R  build.R  project.R
  src/ranges_e4/
    __init__.py   # public exports
    _na.py        # nullable-int helpers (masked arrays)
    rle.py        # Rle (values, lengths); coverage result
    hits.py       # Hits
    iranges.py    # IRanges: construct, accessors, subset, concat, mcols
    seqinfo.py    # Seqinfo + merge
    granges.py    # GRanges: construct, accessors, subset, concat, trim
    _group.py     # split GRanges into (seqlevel[, strand]) groups
    intra.py      # shift narrow resize flank promoters (array kernels)
    inter.py      # range reduce gaps disjoin (array kernels)
    setops.py     # union intersect setdiff + p* versions
    overlaps.py   # overlap kernel + findOverlaps family
    nearest.py    # nearest precede follow distance
    coverage.py
    ordering.py   # order sort rank duplicated unique
  tests/
    conformance/  # Python runner: loader, build, ops, project, compare
    test_conformance.py     # one pytest param per case
    test_case_inventory.py  # clauses <-> cases <-> goldens <-> ADRs
    test_<module>.py        # unit + brute-force oracle tests
```

Group files (one per operation group, same stem in cases/ and goldens/):
`construct`, `intra`, `inter`, `setops`, `overlaps`, `nearest`,
`coverage`, `ordering`, plus `seqinfo`, `trim` for granges.

## 2. Data model (Python)

| Class | Storage | Notes |
|---|---|---|
| `IRanges` | `_start`, `_width`: `np.int64`; `_names`: `np.ndarray[str] \| None`; `_mcols`: `pl.DataFrame` (height = len) | `end = start + width - 1`; width ≥ 0 else `ValueError`; values must fit int32 (R integer) |
| `Seqinfo` | frozen dataclass of tuples: `seqnames`, `seqlengths: int\|None`, `is_circular: bool\|None`, `genome: str\|None` | `merge(a, b)` |
| `GRanges` | `_seqcodes: np.int32` into `seqinfo.seqnames`; `_ranges: IRanges` (holds names, no mcols); `_strand: np.int8` (0 `+`, 1 `-`, 2 `*`); `_mcols`; `_seqinfo` | seqlevels default = order of first appearance (R does this) |
| `Hits` | `query`, `subject`: `np.int64` (0-based); `query_length`, `subject_length` | sorted by query, then subject |
| `Rle` | `values`, `lengths`: numpy | coverage only |

Fixed decisions (implementer writes each as an ADR in stage 2):

- **ADR-0001 indices are 0-based** in the Python API (subscripts, hits,
  order, nearest, revmap). Coordinates stay 1-based closed. Goldens stay
  R-native 1-based; the Python runner adds 1 when projecting index fields.
- **ADR-0002 NA** in integer results (`select="first"` with no hit,
  `nearest` with no subject, `distance` across seqnames) = `np.ma.MaskedArray`
  (mask = NA). Projection maps masked → `null`.
- **ADR-0003 errors** compared by "raised `ValueError`/`TypeError`" only,
  never by message. Warnings recorded in goldens, not compared.
- **ADR-0004 no negative subscripts.** R exclusion (`x[-1]`) is not ported;
  Python negative ints follow numpy. Cases use positive indices only.
- `rank` has keyword-only `ties_method` with **no default** (R errors when
  it is omitted — probed). Not a divergence.

API shape: methods on the classes, snake_case; R argument names map
mechanically `.` → `_` (`min.gapwidth` → `min_gapwidth`,
`ignore.strand` → `ignore_strand`). Op names: `findOverlaps` →
`find_overlaps`, `countOverlaps` → `count_overlaps`, `overlapsAny` →
`overlaps_any`, `subsetByOverlaps` → `subset_by_overlaps`, `isCircular`
→ `is_circular`; others keep their R name (`reduce`, `pintersect`, ...).
Concatenation: `IRanges.concat(*xs)` / `GRanges.concat(*xs)`.

## 3. Core algorithms

Seam: every algorithm is a kernel over plain `start`/`end` int64 arrays
(IRanges logic). GRanges methods split by group key (`_group.py`:
seqlevel × strand, or seqlevel only when `ignore_strand=True`), call the
kernel per group, and reassemble with the group's seqname/strand. No
algorithm is written twice.

| Group | Kernel |
|---|---|
| intra | vectorised arithmetic. GRanges: `-` strand flips `fix`/`start` side for `resize`, `flank`, `promoters`; `shift`/`narrow` ignore strand |
| range/reduce | stable sort by (start, end); sweep merging while `next.start <= run_end + min_gapwidth`; revmap = original indices per run in sorted order; zero-width handling per `drop.empty.ranges` (pin with cases) |
| gaps | complement of reduce within `[start, end]` bounds; GRanges uses seqlengths, and emits whole-sequence gaps for unused strands/seqlevels as R does (pin with cases) |
| disjoin | breakpoints = all starts ∪ all ends+1; keep segments covered ≥1 |
| set ops | `union = reduce(c(x,y))`; `intersect`, `setdiff` via reduce + gaps; p* element-wise (`punion` `fill.gap`, `pintersect` `resolve.empty`) |
| overlaps | sort subject by start, `np.searchsorted` to bound candidates by `start <= q_end + maxgap`, then filter on end/type/minoverlap. Correctness over speed; O(n·m) worst case accepted. GRanges: per seqlevel, then strand compatibility (`*` matches all) unless `ignore_strand` |
| select | first/last = min/max subject per query; arbitrary = any hit (compare mode `oneof`) |
| nearest | precede/follow via sorted starts/ends + searchsorted; nearest = overlap first, else min distance; ties per `select` |
| distance | `max(0, max(starts) - min(ends) - 1)`; zero for overlap/adjacent per R (pin with cases) |
| coverage | diff array (+w at start, −w at end+1), cumsum, run-length encode; `width`, `shift`, `weight`. GRanges → `dict[str, Rle]` over seqlevels, length = seqlength or max end |
| ordering | `np.lexsort` stable; IRanges key (start, width); GRanges (seqcode, strand, start, width); duplicated/unique on the same keys |
| seqinfo merge | union of seqlevels, x first then new ones from y; conflicting non-null length/circularity/genome → `ValueError` (pin R's ordering rules with cases) |
| trim | clamp to `[1, seqlength]` on non-circular seqlevels with known length |

## 4. Case artifacts

### 4.1 Case file (hand-authored, YAML)

```yaml
# cases/iranges/overlaps.yaml
cases:
  - id: ir-fo-zero-width-subject
    clause: IR-OVL-ANY-3
    inputs:
      x: {IRanges: {start: [1, 5], width: [10, 3]}}
      s: {IRanges: {start: [6], width: [0]}}
    call: {op: findOverlaps, args: [x, s]}
  - id: ir-fo-select-arbitrary
    clause: IR-OVL-SELECT-4
    inputs: {x: ..., s: ...}
    call: {op: findOverlaps, args: [x, s], kwargs: {select: arbitrary}}
    compare: {fields: {values: oneof}}
    oneof_from: {kwargs: {select: all}}
    adr: ADR-0005
```

Rules:
- `id`: globally unique, prefix `ir-`/`gr-`. `clause` must exist in `clauses.md`.
- `inputs`: name → one constructor. Allowed: `IRanges{start,end,width,names,mcols}`,
  `GRanges{seqnames,start,end,width,strand,names,mcols,seqinfo}`,
  `Seqinfo{seqnames,seqlengths,isCircular,genome}`, and literal vectors
  `int: [...]`, `lgl: [...]`, `chr: [...]` (e.g. subscripts). `null` = NA.
  `mcols: {col: [...]}`, type inferred from first non-null.
- `call.args`: positional, input names only. `call.kwargs`: primitives,
  R argument names. `{ref: name}` passes an input as a kwarg.
- Pseudo-ops: `identity` (project an input — construction cases), `[`,
  `c`, `length`, and accessors (`start`, `end`, `width`, `names`,
  `seqnames`, `strand`, `seqlevels`, `seqlengths`, `isCircular`, `genome`).
- No expected values in case files. Ever.
- Construction errors (bad width etc.): put the bad spec in `inputs`;
  R records `kind: error`; Python must raise while building.

### 4.2 Golden file (R-generated, JSON)

```json
{"generated_by": {"R": "4.6.0", "IRanges": "2.47.2", "GenomicRanges": "1.65.1", "S4Vectors": "0.51.6"},
 "cases": {"ir-fo-zero-width-subject":
   {"kind": "Hits", "query": [1, 2], "subject": [1, 1], "query_length": 2, "subject_length": 1}}}
```

Projection by result class (same function in both runners):

| R class | kind | fields |
|---|---|---|
| IntegerRanges | `IRanges` | start, end, width, names, mcols{col: vector \| list of int vectors} |
| GRanges | `GRanges` | seqnames (chr), start, end, width, strand (chr), names, mcols, seqinfo{seqnames, seqlengths, isCircular, genome} |
| Hits | `Hits` | query, subject, query_length, subject_length |
| Rle | `Rle` | values, lengths |
| RleList | `RleList` | names, elements[{values, lengths}] |
| Seqinfo | `Seqinfo` | seqnames, seqlengths, isCircular, genome |
| atomic | `vector` | type, values, names |
| condition | `error` | message (informational) |

`oneof` cases also store `candidates`: per element, the allowed values,
computed by R from `oneof_from` (call re-run with overridden kwargs).
Vectors are always JSON arrays (even length 1); NA → `null`.

### 4.3 Compare modes

| Mode | Scope | Meaning | ADR? |
|---|---|---|---|
| `exact` | field | default; same length, values, order, NA positions | no |
| `ignore` | field | not compared (implementation-defined) | yes |
| `oneof` | field | each element ∈ R's candidates | yes |
| `multiset` | rows | row tuples equal ignoring order | yes |
| `set` | per element of list field | e.g. revmap members ignoring order | yes |

Syntax: `compare: {rows: multiset, fields: {names: ignore}}`. Any
non-exact mode without `adr:` fails `test_case_inventory.py`.

Python divergence: case keeps the R golden, adds
`python: {expect_error: true}` or `python: {skip: ...}` plus `adr:`.
Never edit a golden to match Python.

### 4.4 Runners

R (`runners/r/run_cases.R`, one Rscript process for all files):
- `Rscript runners/r/run_cases.R` — write goldens for all case files.
- `... --check` — recompute, diff against committed goldens, exit 1 on
  drift (upstream behavior change detector; maintenance replay entry point).
- `build.R`: spec → object; normalise yaml lists with `NULL` → typed
  vector with `NA`; coerce coordinates `as.integer`. `suppressPackageStartupMessages`.
- Every build + call in `tryCatch`/`withCallingHandlers` → `error`/`warnings`.
- `jsonlite::toJSON(auto_unbox = FALSE, na = "null", digits = NA, pretty = TRUE)`,
  scalars via `unbox`; case order preserved → stable diffs.

Python (`tests/conformance/`):
- `ops.py`: explicit table R op → (callable, index-result flag). Index
  results get +1 on projection; `[` index inputs get −1 when building.
- `project.py` mirrors §4.2. `compare.py` implements §4.3.
- Missing golden → test fails with "run `just goldens`".

`justfile`: `goldens`, `check-goldens`, `test` (`uv run pytest`), `all`.

### 4.5 Inventory checks (`test_case_inventory.py`)

- every clause in `clauses.md` has ≥1 case; every case clause exists;
- every case has a golden; no golden without a case;
- every non-exact compare / `python:` override names an ADR file that exists;
- every TARGET.md operation appears in some clause (list kept in the test).

### 4.6 `clauses.md`

Table `id | statement (one line) | R ref (man page)`. Ids
`IR-<GROUP>-<OP>-<n>`, `GR-...`, `SI-...`. A clause = one behavior worth
stating (e.g. "IR-INTER-REDUCE-2: `min.gapwidth=0` keeps adjacent ranges
separate"). Append-only; retire by marking `(retired)`.

### 4.7 Maintenance

- New behavior: add clause → add case → `just goldens` → implement.
- R upgrade: `just check-goldens`; each drift = re-run `just goldens`,
  commit goldens alone, then fix Python (or ADR).
- Golden diffs are reviewed like code; surprising R output → NOTES.md.

## 5. Work breakdown (commit per line; ≥2 commits per stage)

Each stage: (a) clauses + cases, run R, commit `cases+goldens: <group>`;
(b) implement + unit tests to green, commit `impl: <group>`.
Target 8–15 cases per group; breadth before depth (4 h ceiling).

| # | Stage | Content |
|---|---|---|
| 1 | Harness | `clauses.md`, `justfile`, R runner (build/project/write/check), 3 seed IRanges construct cases, goldens. Verify NA/len-1/empty-vector yaml handling |
| 2 | Python harness + IRanges core | `tests/conformance/*`, inventory test, ADR-0001..0004; IRanges construct/accessors/names/subset/concat/mcols |
| 3 | IRanges intra | shift, narrow, resize (fix start/end/center), flank (start, both), promoters |
| 4 | IRanges inter | range, reduce (min.gapwidth, drop.empty.ranges, with.revmap), gaps (start/end), disjoin |
| 5 | IRanges setops | union, intersect, setdiff, punion, pintersect, psetdiff |
| 6 | IRanges overlaps | findOverlaps types × select, maxgap, minoverlap; count, overlapsAny, subsetByOverlaps; ADR-0005 `oneof` |
| 7 | IRanges nearest | nearest, precede, follow (select first/all/arbitrary), distance |
| 8 | IRanges coverage | Rle; width, shift, weight |
| 9 | IRanges ordering | order, sort (decreasing), rank (ties), duplicated, unique |
| 10 | Seqinfo + GRanges core | Seqinfo, construct, accessors, subset, concat (seqinfo merge + conflict errors), trim |
| 11 | GRanges intra | strand-aware resize/flank/promoters; shift/narrow |
| 12 | GRanges inter + setops | per (seqname, strand); ignore.strand; gaps with seqlengths |
| 13 | GRanges overlaps + nearest | strand compatibility, ignore.strand, distance across seqnames → NA |
| 14 | GRanges coverage + ordering | RleList per seqlevel; order across seqlevels and strand |
| 15 | Edge batch | empty inputs, zero-width everywhere, int32 bounds, circular seqlevels, NA in names/mcols |
| 16 | Close-out | README (design + artifact maintenance), inventory audit, NOTES |

Unit tests beyond cases: brute-force oracles on seeded random inputs for
overlaps, reduce, disjoin, coverage, nearest (pure-Python O(n·m) vs kernel).

## 6. Risks

| Risk | Mitigation |
|---|---|
| R yaml quirks: `[1, null]` → list with NULL; `[]` → `list()`; big ints → double | normalise in `build.R`; stage-1 seed cases exercise each |
| Index base confusion (0- vs 1-based) in runner | single op table flags index results; inventory of index ops reviewed in ADR-0001 |
| Implementation-defined output (arbitrary select, nearest ties, revmap order, result names) | `oneof`/`set`/`ignore` modes, each with ADR — not silent tolerance |
| Projection misses state (e.g. seqinfo dropped, metadata) | project every field in §4.2 for every object result |
| Overfitting to cases | brute-force oracle tests; edge batch (stage 15) |
| Zero-width overlap and `maxgap`/`minoverlap` defaults are subtle (zero-width at 6 hits [1,10] and [5,7]) | dedicated clauses; read `?findOverlaps`; let goldens decide |
| GRanges gaps/seqinfo merge ordering rules | pin with cases before implementing |
| Time ceiling | breadth-first stages; edge batch optional |
| Goldens are R-version specific | `generated_by` header; `--check` drift mode |
