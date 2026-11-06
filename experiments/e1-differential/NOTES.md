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
