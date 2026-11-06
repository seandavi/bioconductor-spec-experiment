# E0 control — notes (append-only)

Judgment calls, deviations from R and human interventions, newest last.
Format: `## YYYY-MM-DD HH:MM — topic`.

## 2026-11-05 22:35 — planner: scaffold and baseline decisions

- Reference stack on this machine: R 4.6.0, IRanges 2.47.2,
  GenomicRanges 1.65.1, S4Vectors 0.51.6.
- Python pinned to 3.13 (`.python-version`; 3.14 is not installed here),
  `requires-python >= 3.12`. Runtime deps numpy + polars (polars for mcols,
  incl. `List(Int64)` revmap column). Dev deps pytest + hypothesis.
- Planned conventions (see PLAN.md §2.4): 0-based index outputs,
  `NO_MATCH = -1` for missing matches and undefined distances, ranks stay
  1-based, Python negative-index semantics, int64 coordinates.
- R's `rank()` on IRanges errors unless `ties.method` is given explicitly
  (checked ad hoc); plan sets the Python default to `"first"`.
- No human intervention so far.
