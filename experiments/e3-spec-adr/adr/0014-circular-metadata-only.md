# ADR-0014 Circular sequences are metadata only

Status: accepted · Clauses: GR.9, GR.18, SEQI.2

## Context
R treats ranges on circular sequences specially: `findOverlaps` and
`coverage` wrap ranges that run past `seqlength` (e.g. `c:30-35` on a
circular length-30 sequence overlaps `c:1-4`; its coverage covers 1–5 and
25–30). TARGET requires `isCircular` accessors, seqinfo merge and `trim`,
not wrap-around semantics.

## Options
1. Implement wrap-around in overlaps, coverage, nearest, set ops.
2. Store/merge the flag, skip it in `trim`, treat levels as linear.

## Decision
Option 2. Differential fuzzing ran with circular levels excluded; the
divergence is known and listed here.

## Consequences
Out-of-bound ranges on circular levels behave as on linear ones. Anyone
needing wrap semantics must add them to the kernels' GRanges blocks.
