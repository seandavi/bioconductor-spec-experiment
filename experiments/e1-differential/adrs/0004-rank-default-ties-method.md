# ADR-0004: `rank()` has a working default `ties_method`

Status: accepted (divergence from R)

## Context
`rank(<IRanges>)` with no `ties.method` fails in R ("'arg' must be of length 1"):
the method passes the unmatched default vector to `match.arg`. Every explicit
`ties.method` works. Case: `rank-default-ties-method`.

## Options
1. Mirror R: require `ties_method`, raise when omitted.
2. Default to `"average"`, base R's documented default.

## Decision
Option 2. R's failure is an accident of the method signature, and base R's
`rank(c(...))` defaults to `"average"`. Explicit values are tested against R.
`"random"` is unsupported (R's RNG stream cannot be reproduced).
