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
