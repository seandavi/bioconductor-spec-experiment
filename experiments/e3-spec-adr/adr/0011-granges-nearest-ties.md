# ADR-0011 GRanges nearest(select="all") keeps all minimal-distance ties

Status: accepted · Clauses: NEAR.2, GR.15

## Context
IRanges `nearest(all)` is exactly the subjects at minimal distance (R fuzz,
2000 random cases, no differences). GRanges `nearest(all)` in R sometimes
returns fewer ties at non-zero distance (a strand/side-dependent subset;
e.g. `+` query `a:6-7`, subjects `a:1-3:+`, `a:10-12:+` → only `a:10-12`).
The rule is not documented and varies with `*` strands.

## Options
1. Reverse-engineer the subset rule.
2. Define as "all minimal-distance ties over comparable subjects".

## Decision
Option 2. `select="arbitrary"` is unaffected (any tie is valid).

## Consequences
Known divergence from R for GRanges `select="all"` with non-zero-distance ties.
