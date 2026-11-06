# NOTES — e1-differential

Append-only. Timestamped judgment calls and human interventions.

- 2026-11-05 23:30 plan: harness cases in YAML (pyyaml dev dep); Python
  converts to JSON for R (jsonlite). R never parses YAML.
- 2026-11-05 23:30 plan: canonical compare form uses R's 1-based indices;
  adapter post-hooks convert Python results. R encoder stays generic.
- 2026-11-05 23:30 plan: R results cached in committed `r_cache/`, keyed by
  case hash + package versions, so `uv run pytest` runs without R and
  version bumps surface as cache diffs.
- 2026-11-05 23:30 plan: R probes — `reduce` keeps a zero-width range
  outside others (`20-19`) and absorbs one inside; zero-width `5-4` overlaps
  `1-10`, `11-10` does not; `gaps(GRanges)` with seqlengths emits full-length
  gaps on empty strands.
- 2026-11-05 23:55 harness: R-side warning handler records but must NOT muffle.
  Muffling changed R's control flow (`IRanges(start=2147483647, width=2)` with
  doubles errors only if the coercion warning propagates).
- 2026-11-05 23:55 R quirk: `IRanges(start=<int vec>)` alone builds runs of
  consecutive +1 integers (and TRUE-runs for logical); implemented, not an ADR.
- 2026-11-05 23:58 retracted: the apparent integer-overflow divergence was a
  harness artifact (warning muffling). R errors for integer inputs too; no ADR.
  Lesson: a divergence is only real after checking the harness did not alter R.
- 2026-11-06 00:20 R quirk (matched, no ADR): `reduce(min.gapwidth=0)` merges
  identical empty ranges but silently drops the repeats from `revmap`
  (e.g. 5/0,5/0,5/3 -> revmap [1],[3]). Modelled as "drop output range equal to previous".
  `gaps` and `disjoin` also ignore empty ranges as coverage; `disjoin` still
  breaks other ranges at an empty range's position (probed).
- 2026-11-06 01:10 overlaps: predicates fitted to R on an exhaustive 24x24 pool
  (types x maxgap x minoverlap), see `_kernels/overlaps.py` docstring. Facts found:
  `within` ignores maxgap<=0; minoverlap only applies when >0; `any` uses the
  two-sided gap test (so empty ranges behave); `any` with maxgap>=0 and minoverlap>0 errors.
- 2026-11-06 01:10 overlaps: hit order within a query is engine-dependent. R builds an
  NCList on the shorter side (ties -> subject). nq>=ns: subjects in (start asc, width desc,
  idx) order; nq<ns: index order. `select="arbitrary"` = first hit in that order
  (subject mode) or last-seen = max index (query mode). Matched exactly; no ADR needed.
  `first`/`last` are min/max subject index in both modes.
- 2026-11-06 01:10 ADR-0003: counts are plain arrays; R's query names stripped by the adapter.
- 2026-11-06 01:10 drop.redundant: modelled as "drop (i,j), i>j, if (j,i) is also a hit";
  matches R on symmetric (any/equal) and asymmetric (within/start) grids.
- 2026-11-06 01:50 nearest/precede/follow: ported from R source of the S4 methods (read via
  `selectMethod`): nearest uses findOverlaps(maxgap=0) first (adjacent counts), then
  precede(first)/follow(last); a tie goes to the follow side (arbitrary) or keeps both (all).
  Self mode (`nearest(x)`) uses findOverlaps(drop.self) whose select="arbitrary" is the max
  subject index (not the engine's first). distance() overflow -> NA (masked) + warning, as R.
