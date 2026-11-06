# ADR-0007 reduce(drop_empty_ranges) drops seeds, not joiners

Status: accepted · Clauses: INTER.2–INTER.4

## Context
R: `reduce({[5,4],[6,8]}, min_gapwidth=2, drop_empty_ranges=True)` → `{[6,8]}`
(no flag: `{[5,8]}`), but `reduce({[1,3],[5,4]}, min_gapwidth=3, drop…=True)`
→ `{[1,4]}`: the zero-width range still extends the group's end.
Naive "filter zero-width first" gives `{[1,3]}`; "filter result" gives `{[5,8]}`.

## Options
1. Pre-filter empties (differs from R above).
2. Post-filter (differs from R above).
3. Scan rule: a zero-width range may not seed a group but may join one.

## Decision
Option 3, stated as INTER.4. Matches both R observations.

## Consequences
Kernel reduce has one extra branch. Tests cover both example pairs.
