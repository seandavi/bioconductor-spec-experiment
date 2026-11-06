# ADR-0008 `within` ignores max_gap ≤ 0 (R quirk replicated)

Status: accepted · Clauses: OVL.4, OVL.5

## Context
For `type="start"|"end"|"equal"`, `max_gap` 0 and -1 both mean "exact". For
`within`, R treats `max_gap ≤ 0` as "no limit" and `max_gap ≥ 1` as a bound
on `width(subject) - width(query)`; `max_gap=0` is therefore *looser* than 1.

## Options
1. Replicate (R is normative).
2. "Fix": 0 = identical widths.

## Decision
Option 1. Callers who want a tight bound pass `max_gap ≥ 1`.

## Consequences
Surprising non-monotonic behavior; documented in OVL.5 with an example.
