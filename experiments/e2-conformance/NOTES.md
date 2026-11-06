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
