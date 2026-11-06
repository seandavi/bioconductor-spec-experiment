# ADR-0004 snake_case names, methods on the first operand

Status: accepted · Clauses: CONV.9

## Context
R generics are `findOverlaps(query, subject)`, args like `min.gapwidth`.

## Options
1. Copy R names (non-PEP8, dots invalid).
2. snake_case; free functions.
3. snake_case; methods on first operand; `concat` module function.

## Decision
Option 3: `query.find_overlaps(subject)`, `x.union(y)`, `x.reduce(min_gapwidth=1)`.
R `c()` becomes `concat(*objs)`. `select`/`type` string values keep R's spelling.

## Consequences
One spelling per concept; mapping is mechanical (CONV.9).
