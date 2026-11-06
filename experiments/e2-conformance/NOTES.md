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
