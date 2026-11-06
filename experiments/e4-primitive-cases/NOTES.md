# NOTES (append-only, timestamped)

- 2026-11-05 planner: scaffolded uv project (`ranges_e4`, numpy + polars; dev: pytest, pyyaml). Wrote PLAN.md.
- 2026-11-05 planner: R probe findings used in PLAN: `rank()` on IRanges errors without explicit `ties.method`; GRanges seqlevels follow order of appearance; zero-width subject at 6 hits [1,10] and [5,7]; `reduce(with.revmap=TRUE)` revmap follows `order(x)`; YAML `null` inside a vector reaches R as `list(..., NULL, ...)`; package-built-under warnings go to stderr.
