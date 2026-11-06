# NOTES (append-only, timestamped)

## 2026-11-05 — planning

- Format: JSON, not YAML. YAML 1.1 coerces `no`/`yes`/`NA`-like tokens; JSON parsers on both sides keep `1` vs `1.0` (jsonlite checked with `simplifyVector = FALSE`).
- Goldens in `expected/`, separate from authored `cases/`. R rewriting authored files would churn diffs.
- Clauses co-located in case files, not a separate registry. Less drift.
- Format indices 1-based (R-native); Python API 0-based; conversion only in `conformance/ops.py`.
- mcols as polars DataFrame (project default dataframe lib; typed nulls; list columns for revmap). Cost: 0-column frames need a `None` sentinel.
- Error cases compare only "raised" vs "not raised", not messages. Warnings recorded, not compared.
- Probed: `findOverlaps` returns `SortedByQueryHits`; `reduce(with.revmap=TRUE)` revmap is `CompressedIntegerList`; `coverage(GRanges)` is `SimpleRleList`; `mcols(IRanges)` is NULL when empty; Seqinfo is its own package (1.3.0); zero-width IRanges at 5 overlaps [2,11] under defaults.
- R prints "built under R 4.6.1" warnings on attach; runners must suppress.

## 2026-11-05 — Phase A (suite)

- serialization: FORMAT.md + R codec + custom JSON writer + selftest took ~25 min of the session, no redesign. jsonlite's `toJSON` not used for output: it can't keep `NA` vs `NaN` apart in lists and pretty-prints every array element on its own line.
- Dropped plan item: "mixed shorthand" is an error rather than coerced.
- R surprises pinned in iranges-core (R is right, Python must copy):
  - `IRanges(a, b)` positional 2nd arg is `end`.
  - Fractional start/width are truncated, not rejected (`IRanges(1.5, 2)` → start 1).
  - `names` shorter than ranges is NA-padded; longer is an error.
  - `x[integer(0)]` keeps `names = character(0)` (not NULL).
  - `c()` of named + unnamed pads names with `""`; of mcols + no-mcols pads columns with NA.
- findOverlaps: Hits for one query come out in NCList traversal order (e.g. `3-7 3-3 3-4`), not ascending subject index. Not a contract → compare mode `hits-query-grouped` (queries exact, subjects per query as a set). `select=first/last` is index-based and exact.
- nearest(select="arbitrary") breaks ties unspecifiably; fixtures avoid ties for it. `follow` takes select last|all, `precede` first|all (opposite pick errors).
- p* set ops do not recycle length-1 (error); `distance` does (warns on fractional). `psetdiff` errors whenever y lies strictly inside x.
- coverage: positions < 1 silently dropped; NA weight gives garbage (not pinned).
- reduce/disjoin `revmap` order = indices sorted by (start, end), ties by index.

## 2026-11-05 — Phase A done (suite: 884 cases, 16 groups; `run.R check` green)

- Time: suite authoring ~2.5 h wall incl. reading goldens; codec/runner ~25 min of that. Reading R goldens for surprises was the dominant cost, and it paid: most quirks below were not guessable.
- Case counts: iranges core 66 / intra 69 / inter 57 / setops 50 / overlaps 133 / nearest 42 / coverage 29 / order 31; seqinfo 45; granges core 80 / intra 57 / inter 38 / setops 48 / overlaps 62 / nearest 28 / coverage 21 / order 26.
- More R quirks pinned (R right, Python must copy):
  - `reduce` merges identical zero-width ranges even with `min.gapwidth=0`; other zero-width ranges follow `next.start - max.end.excl < min.gapwidth` (fuzzed vs R, 0 mismatches over 4000 draws).
  - `gaps(x)` default window = [min start, max end] of the *non-empty* reduced ranges.
  - `findOverlaps(type="within", maxgap<=0)` is unconstrained; maxgap>0 bounds `width(s) - width(q)`. `type="start"/"end"/"equal"` need no overlap at all, only |start diff| <= max(maxgap,0).
  - `GRanges precede/follow`: for x `*` the direction is taken from the *subject* strand (`-` subject reverses).
  - `pintersect(GRanges)` adds mcols `hit`; strand = the non-`*` of the pair; `punion` keeps x's strand; incompatible pairs error (`punion`) / flag (`pintersect`) / pass x through (`psetdiff`).
  - `rank(IRanges)` without `ties.method` errors in R (unevaluated default vector); `precede/follow(GRanges, empty subject)` return Hits not integer. Both look like R bugs → candidates for xfail.
  - `GRanges(names=...)` is an mcols column; `GRanges(mcols=DFrame)` makes columns `mcols.<col>` (case dropped).
  - `order(GRanges, ignore.strand=TRUE)` errors; `sort`/`rank` accept it.
- Compare mode `hits-query-grouped` declared at file level (`"compare"` key) for overlaps/nearest groups.

## 2026-11-05/06 — Phases B–D (Python harness, implementation, close-out)

Timing correction: the Phase A entry above gave wall-time *estimates* ("~2.5 h"). Git commit
timestamps (the only clock I trust) say: scaffold→suite done 23:27→23:47 (20 min), Python
core+codec+adapters+runner 23:47→23:57 (10 min), tests/types/report/README ~15 min more. Wall time
is not comparable to the estimates; use commits and counts, not those numbers.

### Result

- Suite: 882 cases / 145 clauses / 17 groups. R `check` 882/882.
- Python: 878 pass, 4 xfail, 0 fail on the first complete run after two kernel fixes (below).
  `uv run pytest`: 1951 passed, 4 xfailed. `mypy --strict src` clean.
- Mutation spot-checks (flip resize-center rounding, drop the zero-width reduce rule, loosen
  `within` maxgap, ignore strand-merge in p*): the suite fails 9 / 1 / 22 / 4 cases. It does catch regressions.

### What the first Python run found (all four were my misreadings of R, not suite bugs)

1. `disjoin`: zero-width ranges add no coverage but DO split the range they sit inside (`1/10` + `3/0` → `1/2, 3/8`).
2. `pintersect(GRanges)`: compatible pairs that do not overlap are `hit=FALSE` (zero-width at x's start), not an error.
3. Two nearest cases used `select="arbitrary"` on tied data → the suite pinned an unspecified choice. Removed the cases (`ir.nearest.reversed`, `gr.nearest.reversed`) and moved `ir.nearest.names-dropped` to a tie-free subject. Lesson: the generator happily pins anything; clause review must ask "is this a contract?".
4. Quirks that R gets wrong or inconsistent → `xfail.json` (strict): `ir.rank.default-is-error`, `gr.precede.empty-subject`, `gr.follow.empty-subject`, `gr.nearest.empty-subject`.

### Deviations from PLAN.md

| Plan | Done | Why |
|---|---|---|
| compare mode `hits-unordered` | `hits-query-grouped`, settable per file | R emits sorted queries but NCList-order subjects; pinning only what is contractual |
| decode inputs straight to Python types | `Vec(kind, values, names)` for atomic args | JSON → Python loses int/double, NA, names; adapters need R's view |
| `RleList` class | `dict[str, Rle]` | one use (GRanges coverage); a class adds nothing |
| gate `uv run pytest` with a marker | no gating; suite green from the first commit that wires it | there was never a long red phase (core written before the runner) |
| order of work 7–20 (op groups, one commit each) | core + kernels written in one pass, then the runner | the suite already existed, so a red/green loop per group added nothing |
| nearest via sorted search | per-query numpy scan | simpler, obviously right; documented as O(nq·ns) |

### serialization: where it was hard and what it cost

| Item | Cost | Notes |
|---|---|---|
| Format design + R codec + JSON writer + selftest | ~1 commit, no redesign | custom JSON writer instead of `toJSON`: it cannot keep NA/NaN apart in lists and prints one array element per line |
| Python codec | ~1 commit; **797/797 goldens round-tripped on the first run** | the tagged format left no ambiguity; `Vec` was the one addition |
| Per-op adapters (`ops.py`, ~750 lines, largest file of the harness) | biggest harness cost | R conventions leak here: 1-based/negative/zero subscripts, `NA` args (`narrow(start=NA)`), `ignore.strand` only for GRanges, `...` varargs (`range`, `order`), named vectors (`seqlengths<-`, `coverage(width=)`), partial R defaults |
| NA in numeric inputs | design item | Python ints have no NA → `Vec` keeps `None`; outputs use `NA_INDEX=-1` and convert at the boundary |
| Empty vectors | design item | `[]` has no type → tags are mandatory in inputs for empties |
| mcols | polars cannot hold 0 columns with rows → `None` sentinel + `nrow` carried by the owner; recycling rule (cyclic, even non-divisors) is R's |
| `revmap` | 0-based `List(Int32)` column in Python vs 1-based `IntegerList` in the format; the codec special-cases the column name `revmap` | hacky but contained |
| `Rle` of `factor` | `seqnames()`/`strand()` return `Rle(Factor)`; needed `Factor` + `Rle.from_factor` | |
| `Hits` | order inside a query is an R artifact → compare mode | cost: one design iteration and a NOTES entry |
| Named results | `seqlengths(x)` etc. are named vectors → names are part of the tree | |
| Warnings / error messages | not compared | cheap and safe; means the suite cannot see "wrong warning" regressions |

Biggest unexpected cost was not serialization but **reading goldens to find which R behaviours are
contracts and which are accidents** (~25 quirks listed in the Phase A entries).

### Judgment calls

- `select="arbitrary"` → lowest index. `nearest` tie-breaking otherwise unspecified.
- `GRanges.coverage` default length for NA seqlength = max shifted end of that level's ranges (0 if none).
- `GRanges.pintersect` for compat+non-overlap: strand of the pair, `hit=FALSE`; for incompatible: x's strand. Only the second is pinned by a case; the first follows from one R probe.
- `seqlevels<-` on data: rename if nothing in common and same length, else map by name (drop unused / add new). Matches every case; other partial-overlap forms are untested.
- R locale: levels are first-appearance order, so no case depends on collation.
- Python `IRanges.rank()` defaults to `"average"` (R's own default errors).
