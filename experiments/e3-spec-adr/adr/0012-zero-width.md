# ADR-0012 Zero-width ranges follow R's position-between-bases model

Status: accepted · Clauses: OVL.3, INTER.3, INTER.7, SET.3–SET.5, COV.2

## Context
`[s,s-1]` sits between positions `s-1` and `s`. R treats it as: overlapping
a range only if strictly inside; splitting a range in `disjoin` and
`psetdiff`; covering nothing in `coverage`/`gaps`/set operations; still
counting for `range`, `reduce` seeds, and coverage length.

## Options
1. Treat zero-width as empty/ignored everywhere (simpler; wrong vs R).
2. Replicate each operation's R behavior, specced per clause.

## Decision
Option 2. No global rule exists, so each clause lists its zero-width behavior
and examples.

## Consequences
Each kernel has explicit zero-width tests.
