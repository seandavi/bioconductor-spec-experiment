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
