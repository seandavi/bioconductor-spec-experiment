# NOTES — e3-spec-adr

Append-only. Timestamped judgment calls.

- 2026-11-05 planner: element indices 0-based, coordinates 1-based (ADR-0001 to record).
- 2026-11-05 planner: NA integer results as numpy masked arrays; keeps core numpy-only, polars only for mcols.
- 2026-11-05 planner: traceability test (every clause cited by ≥1 test) added to keep spec and tests linked; not a conformance suite.
- 2026-11-05 planner: no R goldens committed; spec examples are the record of R behavior.
- 2026-11-05 implementer: PLAN assumed S1–S4 spec stages were done; spec/ and adr/ were empty, so I ran S1–S4 first (R-probed every clause), committed spec + ADR-0001..0012, then I1..I6.
- 2026-11-05 implementer: added a gitignored scratch differential fuzz (spec/scratch/fuzz_gen.R + fuzz_cmp.py; R emits cases+answers, Python recomputes). Not a golden/harness artifact; it only guided spec amendments. It found 5 corner cases (ADR-0013, INTER.3/5/6, INTRA.2, NEAR.2); every amendment is in the commit that fixed it.
- 2026-11-05 implementer: deviation from PLAN algorithms: overlaps use a candidate window from sorted subject starts + max subject width (chunked), nearest is built on overlaps(max_gap=0) plus two searchsorted groups; reduce is a vectorised cummax scan with a Python-loop fallback only for drop_empty_ranges.
- 2026-11-05 implementer: `nearest(subject=None)` for zero-width x: R returns self (twice); left implementation-defined in NEAR.2 rather than replicated.
- 2026-11-05 implementer: GRanges `nearest(select="all")` tie subset not replicated (ADR-0011).
- 2026-11-05 implementer: GRanges fuzz found R behaviours I would not have guessed (all now in spec): intersect/setdiff clip to seqlength but union does not (GR.12); parallel-op strand rules differ per op (GR.13); `*` queries merge `+/*` and mirrored `-` neighbor candidates by distance (GR.15); GRanges promoters has no ignore.strand; GRanges pintersect has no resolve.empty.
- 2026-11-05 implementer: decided circular wrap-around is out of scope (ADR-0014): TARGET lists only isCircular/merge/trim. R wraps in findOverlaps/coverage; fuzz ran with circular=FALSE.
- 2026-11-05 implementer: not replicated and not specced: R GRanges precede/follow with a zero-width subject at position 1 (offset quirk); self-nearest of zero-width x. Both marked implementation-defined in the spec.
- 2026-11-05 implementer: closed out — all TARGET operations public on IRanges and GRanges; trace allow-list empty; ~110k scratch R-vs-Python comparisons clean apart from ADR-0011 (GRanges nearest all) and the quirks above.
