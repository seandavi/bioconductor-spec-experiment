# ADR-0003: atomic results are plain arrays without names

Status: accepted (divergence from R)

## Context
R returns some vectors named after the query (`countOverlaps(x, y)` carries
`names(x)`). Python results are numpy arrays; names live on the object (`x.names`).

## Options
1. Return a named wrapper (dict, or ndarray subclass) from every vector-valued op.
2. Plain arrays; the user zips with `x.names` when needed.

## Decision
Option 2. A wrapper type would leak into every downstream numpy use for a
convenience R users get mostly by accident. The adapter strips R names for the
functions listed in `adapter.R_UNNAMED` (currently `countOverlaps`); cases still
assert the values. Names on object results (IRanges/GRanges) are kept.
