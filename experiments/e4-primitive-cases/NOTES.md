# NOTES (append-only, timestamped)

- 2026-11-05 planner: scaffolded uv project (`ranges_e4`, numpy + polars; dev: pytest, pyyaml). Wrote PLAN.md.
- 2026-11-05 planner: R probe findings used in PLAN: `rank()` on IRanges errors without explicit `ties.method`; GRanges seqlevels follow order of appearance; zero-width subject at 6 hits [1,10] and [5,7]; `reduce(with.revmap=TRUE)` revmap follows `order(x)`; YAML `null` inside a vector reaches R as `list(..., NULL, ...)`; package-built-under warnings go to stderr.
- 2026-11-05 impl: stage 1 done. R runner is generic: `do.call(call$op, args+kwargs)` so R op names are used verbatim (no op table in R). YAML gotcha confirmed: R yaml parses unquoted `Y`/`n`/`yes`/`on` as logicals, PyYAML does not (only yes/no/on/off/true/false). Rule: quote any string that could be a YAML 1.1 boolean; use `chrX` style seqnames. Ints > int32 become NA in R yaml: keep coordinates within int32.
- 2026-11-05 impl: R-yaml `y`/`n` bool quirk bit immediately (input named `y` became TRUE). Fixed in run_cases.R with yaml handlers keeping y/Y/n/N as strings; `yes/no/on/off/true/false` stay booleans in both parsers. No lint needed.
- 2026-11-05 impl: stage 2 done. test_case_inventory target-ops check is vacuous until a STAGE_DONE file exists; remove that escape hatch in stage 16 (close-out).
- 2026-11-05 impl: R probes pinned: gaps() defaults use non-empty coverage extent, zero-width ranges never split gaps; reduce zero-width rule = same sweep rule as non-empty (merge iff start-runend-1 < min.gapwidth); disjoin ignores zero-width. Oracle tests in tests/test_inter.py.
- 2026-11-05 impl: plan deviations (all committed with the stage that made them):
  - Stage 3: kernels take a `reverse` mask from the start, so GRanges strand-awareness (stage 11) reused them unchanged.
  - Stage 6: overlap kernel is a per-query sorted-prefix filter, not an interval tree; GRanges reuse it via a 2^34 per-seqlevel coordinate offset and 7 strand-combination calls.
  - Stage 10: `si-` case ids became `gr-si-` (inventory test requires ir-/gr- prefixes). R projects `seqnames()`/`strand()` as chr vectors (ADR-0008).
  - Stage 12: GRanges `punion/pintersect/psetdiff` were not in the plan; added because TARGET lists `p*` ops. R rules reverse-engineered (hit mcol, strand resolution) and pinned with cases.
  - Stage 13: precede/follow direction rule for `*` queries (take the subject's strand) found by probing; the plan assumed plain per-group kernels.
  - Stage 15: edge probing found R defects at end=2^31-1 (gaps, setdiff, disjoin, maxgap) and a Hits return for zero-length precede/follow: skipped cases + ADR-0010; Python pinned in tests/test_edge.py.
- 2026-11-05 impl: R behaviours worth knowing (all pinned by cases): `findOverlaps` within + `maxgap=0` means no slack limit (only `maxgap>0` limits); `gaps()` on GRanges visits every seqlevel x strand with start=1,end=seqlength; `coverage()` wraps ranges on circular sequences; psetdiff on GRanges drops mcols while pintersect keeps them and adds `hit`; `*` is not a wildcard in sequence set ops but is in overlaps/nearest.
- 2026-11-05 impl: R probes were run ad hoc in /tmp scripts (not committed) and only the resulting behaviours were turned into clauses + cases. Recommendation for a next iteration: keep a `probes/` dir of those scripts so clause provenance is reproducible.
- 2026-11-05 impl: stage 16 close-out: removed STAGE_DONE escape hatch from the inventory test; TARGET operation coverage is now strict, and checked per class (cases) and per API (test_public_api.py).
