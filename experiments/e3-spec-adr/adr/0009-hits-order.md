# ADR-0009 Hits sorted by (query, subject)

Status: accepted · Clauses: OVL.7, NEAR.2

## Context
R returns hits grouped by query, but subjects within a query follow the
interval-index traversal (start ascending, end descending), e.g. subjects
`0-4-2-1-3` for one query. That order is an implementation detail of R.

## Options
1. Reproduce R's traversal order.
2. Sort by `(query, subject)`.

## Decision
Option 2: deterministic and cheap to state. Set equality with R holds.

## Consequences
Order-sensitive comparison with R needs sorting on the R side.
