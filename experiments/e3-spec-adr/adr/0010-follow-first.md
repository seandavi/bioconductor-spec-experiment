# ADR-0010 follow(select="first") returns the highest index of a tie

Status: accepted · Clauses: NEAR.3, NEAR.4, GR.15

## Context
When several subjects tie for nearest-following (same largest end), R's
`follow(select="first")` returns the highest index, while `precede(first)`
returns the lowest. `select="all"` lists every tie.

## Options
1. Lowest index for both (symmetric, differs from R).
2. Replicate R.

## Decision
Option 2. Attached to the operation, not to geometry: on GRanges, a `-` query
still gets lowest for `precede`, highest for `follow`.

## Consequences
Asymmetric; spec and tests give explicit ties.
