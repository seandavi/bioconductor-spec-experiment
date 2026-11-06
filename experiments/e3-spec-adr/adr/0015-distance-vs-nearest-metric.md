# ADR-0015 `distance` and `nearest` use different metrics

Status: accepted · Clauses: NEAR.1, NEAR.2, GR.15

## Context
IRanges 2.47.3 / GenomicRanges 1.65.2 changed `distance()`: overlap gives 0,
otherwise gap width + 1 (adjacent = 1). The release note says `nearest()`,
`precede()` and `follow()` are unchanged. Until then the port used one kernel
(gap width) for both `distance` and the tie-breaking in `nearest`. R's new
`distance` is a wrapper (`ifelse(overlap, 0, old + 1)`) and zero-width ranges
never count as overlapping, so `distance([5,4],[5,4])` is 1.

## Options
1. Change the shared kernel; `nearest` ties between adjacent and overlapping
   subjects would break (R still treats them as tied at 0).
2. Keep a separate gap-width kernel for `nearest` (and GRanges `*` merging) and
   give `distance` its own.

## Decision
Option 2: `_kernels.gap_width` serves `nearest`/GR.15 star merging;
`_kernels.distance` implements NEAR.1. Verified against R: 50k random
IRanges pairs and 50k GRanges pairs (incl. strand/seqlevel NA) match; R
`nearest(select="all")` still ties adjacent with overlapping subjects.

## Consequences
Two near-identical kernels; the difference is the point and is named in the
spec. Supersedes the "adjacent = 0" examples of the earlier NEAR.1 text (the
clause is amended in place).
