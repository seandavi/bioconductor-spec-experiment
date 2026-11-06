# NOTES (append-only, timestamped)

Judgment calls and human interventions. Tag R-behaviour guesses `UNVERIFIED`.

- 2026-11-05 23:25 [plan] Scaffolded uv project (numpy, polars; dev: pytest, hypothesis). Python ≥3.13.
- 2026-11-05 23:40 [plan] PLAN.md written. Key calls: 0-based positional indices, `-1` for R `NA` index, int64 coords, polars mcols, grouped numpy kernels shared by IRanges/GRanges, methods not free functions. R edge-case guesses tagged `UNVERIFIED` in PLAN §4.
- 2026-11-05 23:40 [impl] Implementation stages 1-16 done, one+ commit each. Deviations from PLAN and judgment calls below. All R guesses are `UNVERIFIED` (no R access).
- [impl] PLAN deviation: one name cannot be both classmethod and method, so `x.concat(*others)` + `IRanges.combine(list)` / `GRanges.combine(list)`.
- [impl] PLAN deviation: parallel ops (`punion`, `pintersect`, `psetdiff`) keep x's names and mcols (PLAN said drop). Guess: R keeps x's metadata for pairwise results. `UNVERIFIED`.
- [impl] PLAN deviation: `union/intersect/setdiff` use one sweep over merged breakpoints (`setops.sweep_kernel`) instead of the `gaps(union(gaps(x),y))` formulas. Same normalised output (reduced, non-empty); no lo/hi bookkeeping. Property-tested against position-set brute force.
- [impl] PLAN deviation: `type="within"` with `maxgap >= 0` requires `(qs-ss)+(se-qe) <= maxgap` (my recollection of R's meaning: max slack); `maxgap=-1` means unrestricted. PLAN said maxgap ignored. `UNVERIFIED`.
- [impl] Overlap predicate: empty range at p vs non-empty [a,b] hits (maxgap=-1) iff a < p <= b; touching at a or b+1 only with maxgap >= 0. Two empties never hit at maxgap=-1. `UNVERIFIED`.
- [impl] `minoverlap` is only applied when > 0, so disjoint pairs admitted by a large maxgap with type start/end are not dropped by the default `minoverlap=0`.
- [impl] `narrow(width=k)` alone keeps the leading k positions; `narrow(start=len+1)` gives an empty range at end+1 (no error); error only beyond that. `UNVERIFIED`.
- [impl] `resize(fix="center")` start = `s + (w - w') // 2` (floor). `UNVERIFIED`.
- [impl] `reduce(with_revmap=True)`: indices ascending within each output range (0-based). Order inside a revmap entry in R unknown. `UNVERIFIED`.
- [impl] `pintersect` default errors only when the result width would be negative (gap > 0); adjacent ranges give an empty range. `UNVERIFIED`.
- [impl] `disjoin`: empty ranges add breakpoints (start, end+1) but never appear. `UNVERIFIED`.
- [impl] `gaps` (GRanges): every seqlevel x {+,-,*}; start=1, end=seqlength else max end on the seqlevel, seqlevel skipped if neither. `UNVERIFIED`.
- [impl] Strand rules (GRanges): overlaps/nearest compatible iff equal or either `*`. `type=start/end` flips to the opposite end when either side is `-`. precede/follow: `-` query looks to lower coordinates, `*` as `+`; `ignore_strand` treats all as `+`. Parallel ops: strands equal or one `*` (result takes the non-`*`), else ValueError; `ignore_strand` -> result `*` (psetdiff keeps x's). `UNVERIFIED` throughout.
- [impl] `distance` on GRanges returns -1 (`NO_HIT`) for different seqlevels or incompatible strands (R: NA).
- [impl] `rank` returns 1-based values (they are ranks, not positional indices); float64 for "average". Positional indices elsewhere are 0-based.
- [impl] `select="arbitrary"` is deterministic (= first / lowest index).
- [impl] `findOverlaps` without subject uses self including self-hits; `nearest/precede/follow` without subject exclude self-matches (this is how I recall R).
- [impl] Bug found by hypothesis: coverage default width went negative when all ranges ended before 1. Fixed (own commit).
- [impl] Perf smoke (500k ranges): reduce 0.05s, overlaps 0.7s, nearest 2.1s. Worst case remains one huge subject (candidate window), see README.
- [impl] Not done / not attempted: circular-sequence wraparound; `disjoin` revmap; camelCase R-name aliases; free-function API.
- 2026-11-06 [review] Independent review of HEAD 7322367: 4 findings, all accepted, none rejected (one commit each, regression tests in tests/test_review_fixes.py).
  - F1 (major) GRanges `resize/flank/promoters` lacked `ignore_strand`: added; `ignore_strand=True` treats every range as `+`. Signature test added for all strand-aware GRanges ops.
  - F2 (major) mcols aliasing: DataFrames are cloned on ingestion and returned as clones from `.mcols` (polars clone is cheap, buffers shared copy-on-write). Internal sharing between derived objects is safe because nothing mutates in place.
  - F3 (major) `with_seqinfo` could not drop unused levels: `_remap` now maps only codes in use; removing a used level still raises ValueError.
  - F4 (minor) empty boolean mask of wrong length now raises like any other wrong-length mask; `x[[]]` stays a valid empty subset.
- 2026-11-06 [upstream] IRanges 2.47.3 / GenomicRanges 1.65.2 release note: `distance()` now returns 0 for overlapping ranges and gap width + 1 otherwise (adjacent = 1L); GenomicRanges delegates to IRanges; `nearest/precede/follow` unchanged. Process (README "Tests and how to maintain them"): fix code, fix tests, append here. No R access, so this is taken from the release note only, not run in R (`UNVERIFIED` against R). Changed `nearest.distance` (`np.where(gap < 0, 0, gap + 1)`); updated example tests (IRanges `[3,0,0,0]` -> `[4,0,0,1]`; GRanges `[4,0,NA,4]` -> `[5,0,NA,5]`, ignore_strand `4` -> `5`), added an adjacency/symmetry test, README and PLAN §Nearest note. Judgment: empty ranges still follow the old overlap rule (`gap < 0` is overlap); the note does not say otherwise. `nearest()` keeps ranking by gap width (internal, unchanged), so adjacent subjects still win ties as before.
