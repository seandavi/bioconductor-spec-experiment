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

## 2026-11-05 23:20 — implementer: stages 1–5 findings and deviations from PLAN.md

- Constructor: `IRanges(5)` (start only) = width 1 (R behaviour); inconsistent
  start/end/width errors; only length-1 args are recycled (R also recycles
  other lengths with a warning — not reproduced). Fractional floats rejected.
- Zero-column polars frames have no height, so "no mcols" is `pl.DataFrame()`
  and `take_rows` is a no-op for it.
- `narrow`: offsets are relative to each range and never widen it (R's
  `allow.nonnarrowing=FALSE`); `width` needs exactly one of start/end.
- `resize(fix="center")`: `start + floor((w - new_w)/2)` (checked odd/even/zero-width).
- **PLAN.md was wrong** about `reduce(with_revmap=True)`: revmap entries come
  in (start, end, index) sweep order, not ascending index (R: `IRanges(c(5,1,3),c(6,2,8))`
  -> revmap `2 3 1`). `disjoin` revmap *is* ascending.
- `reduce`: default `drop_empty_ranges=FALSE`; empty ranges join the sweep, so
  an isolated one survives and one touching/inside a real range is absorbed.
  Merge rule `start <= running_max_end + min_gapwidth`.
- `gaps`: missing bounds default to the first start / last end of the
  *non-empty* reduced ranges; only-empty input with a missing bound -> no gaps.
- `disjoin`: empty ranges contribute a breakpoint at their start (R splits
  `1-4` at an empty `3-2`) but are never covered.
- Set ops (`union/intersect/setdiff`) are position-set algebra returning
  reduced non-empty ranges; no names/mcols. `p*` ops copy R's rules (checked
  from the R source for the formulae, adjacent pairs give a definite empty
  `pintersect`).
- `rank` default `ties_method="first"` (R: must be specified); `order(decreasing=True)`
  keeps ties in original order, like R.

## 2026-11-05 23:55 — implementer: overlap predicate (derived empirically in R)

The IRanges C code isn't readable here, so the `findOverlaps` predicate was
inferred from small ad hoc R tables (one subject, sweep of query starts/widths)
and is documented in `_overlaps.py`. Notable points:
- `type="any"`: closed-interval overlap `qs <= se+g+1 and ss <= qe+g+1`; this
  makes a width-0 query at p overlap `[a,b]` iff `a < p <= b`.
- `maxgap` for start/end/equal = max allowed |difference| (-1 acts as 0).
- `type="within"`: `maxgap` bounds how much wider the subject may be, but
  **maxgap=0 acts as unlimited** (same as -1) — R quirk, reproduced and tested.
- `minoverlap` is an extra conjunct (`overlap length >= m`) for every type.
- Hits order: R does *not* sort subjects within a query (e.g. self-overlaps
  give 1,3,2,5 for query 1). We always sort by (query, subject) and treat that
  as a deliberate deviation; tests compare as sorted pairs.
- `select="arbitrary"` returns the first hit.
- Hypothesis brute-force test caught a kernel bug (candidate bound for
  `type="end"` with empty subjects); fixed.
