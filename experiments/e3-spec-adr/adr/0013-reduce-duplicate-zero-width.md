# ADR-0013 reduce(min_gapwidth=0) drops repeated zero-width ranges

Status: accepted · Clauses: INTER.3, INTER.5

## Context
`reduce({[14,13],[14,13],[14,20]}, min_gapwidth=0)` → `{[14,13],[14,20]}`,
revmap `((0),(2))`: R silently loses the second identical zero-width range,
also from `revmap`. With `min_gapwidth ≥ 1` such ranges merge and all are
listed. Found by differential fuzzing against R, not from docs.

## Options
1. Keep every input index in some revmap list (principled, differs from R).
2. Replicate R.

## Decision
Option 2 (R is normative): an unabsorbed zero-width range equal to its scan
predecessor is dropped when `min_gapwidth=0`.

## Consequences
`revmap` is not a partition of the input in this corner. Tested in INTER.3.
