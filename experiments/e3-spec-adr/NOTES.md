# NOTES — e3-spec-adr

Append-only. Timestamped judgment calls.

- 2026-11-05 planner: element indices 0-based, coordinates 1-based (ADR-0001 to record).
- 2026-11-05 planner: NA integer results as numpy masked arrays; keeps core numpy-only, polars only for mcols.
- 2026-11-05 planner: traceability test (every clause cited by ≥1 test) added to keep spec and tests linked; not a conformance suite.
- 2026-11-05 planner: no R goldens committed; spec examples are the record of R behavior.
