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

## 2026-11-06 00:40 — implementer: stages 7–9 (+11 folded in): nearest, coverage, Seqinfo, GRanges core

- IRanges `nearest` follows the R source exactly (maxgap=0 overlap first, else
  closer of precede/follow, ties go to the follow side; `select="all"` keeps
  both). GRanges `nearest` uses a *different* R rule on ties (lower subject
  index wins); to be reproduced in the GRanges stage.
- `distance` is the gap width (0 for overlapping *or adjacent*): read from the R
  source, not the docs prose.
- `coverage`: default length = max end (an only-empty input still gives `end`
  zeros), ranges clipped to [1, length]; implemented as a sparse event sweep.
- GRanges storage deviates from PLAN.md: strand codes are `0=+, 1=-, 2=*` (R's
  level order) rather than +1/-1/0, so sorting by code sorts like R.
- Seqlevels default to *first appearance* (R 4.6/GenomicRanges 1.65 behaviour;
  checked with `c("chr2","chr10","chr1")`), not sorted order.
- GRanges ctor recycles length-1 seqnames/ranges/strand like R (length-1 `ranges`
  with several seqnames is valid).
- Out-of-bound ranges (known, non-circular length) warn with `OutOfBoundWarning`
  at construction/concat/shift/resize/flank/promoters/narrow, not on subsetting.
  Message text differs from R's.
- `trim`: start clipped into [1, len+1], end into [0, len]: (150,160) on len 100
  -> (101,100); (-5,-2) -> (1,0). Circular / unknown-length sequences untouched.
- `Seqinfo.merge` warns only when *both* sides have levels missing from the other
  (R wording differs).
- `with_seqinfo` errors if the new seqinfo drops a level in use (R's `seqinfo<-`
  has more modes: pruning, renaming). Not reproduced.

## 2026-11-06 01:30 — implementer: GRanges overlaps / nearest / coverage (stage 13)

Read the GenomicRanges source (`.GenomicRanges_findPrecedeFollow`, `.nearest`) for
the strand logic and confirmed with ad hoc runs:
- `findOverlaps`: same seqname, strands compatible (`*` matches anything) unless
  `ignore_strand`; the IRanges predicate is applied per pair, so maxgap/minoverlap/
  type behave exactly as for IRanges.
- `precede`/`follow`: `+` and `*`-subject look in the + direction for `+` queries;
  `-` queries look the opposite way among `-`/`*` subjects; a `*` query tries both
  (`+`/`*` subjects one way, only `-` subjects the other) and keeps the closer
  side(s) — equal distance keeps both. `select="first"` = lowest subject index
  among ties in either direction, `"last"` = highest.
- GRanges `nearest` differs from IRanges `nearest` on equidistant neighbours:
  `arbitrary` picks the lower subject index; `all` returns only the precede side.
- `distance` returns -1 (NO_MATCH) for different seqnames or `+` vs `-`.
- `coverage` returns a plain dict `{seqlevel: Rle}` (R: RleList) over every
  seqlevel, ignoring strand.
- Hits order: always sorted by (query, subject) — see the IRanges note.
- Not implemented: `distanceToNearest`, `selectNearest`, `overlapsRanges`,
  `findOverlaps` with integer vectors as query, `circle.length`/circular wrap,
  `nthread`.

## 2026-11-06 02:30 — implementer: wrap-up (stage 14)

- All operations in TARGET.md's required scope are public API (`tests/test_scope.py`
  asserts the surface). `uv run pytest`: 300 passed; `mypy src` clean.
- Redesign worth knowing: the first overlap engine bounded candidate subjects only
  from above (by start), which made `type=any/end/within` quadratic (200k x 200k
  never finished). Added a lower bound from the running max of subject ends;
  now ~1 s. The hypothesis brute-force test guards the pruning's correctness.
- Property tests corrected twice by their own failures (my invariants were wrong,
  the code was right): reduce with `min_gapwidth > 1` legitimately covers more
  positions than its input.
- Known unverified spots (expectations written from my reading of R, not each
  checked in R): `Hits`-level behaviour of `select="arbitrary"` (we return the
  first), GRanges `find_overlaps` on seqlevel sets that only partly overlap, the
  exact wording/trigger of R's warnings, `IRanges.nearest` tie rule beyond the
  cases listed above, and any combination of `*`-strand subjects/queries in
  `precede`/`follow` beyond the examples tested.
- No human intervention at any point.
