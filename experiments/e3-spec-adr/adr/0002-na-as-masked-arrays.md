# ADR-0002 Missing integer results are masked arrays

Status: accepted · Clauses: CONV.4, OVL.8, NEAR.2–4, GR.15

## Context
`select="first"`, `nearest`, `precede`, `follow`, GRanges `distance` can be
NA. Plain int64 arrays have no NA.

## Options
1. `float64` with `nan` (loses integer dtype, precision >2^53).
2. Sentinel `-1` (collides with negative subscripts, ADR-0001).
3. `numpy.ma.MaskedArray[int64]`.
4. polars Series with nulls.

## Decision
Option 3: keeps the core numpy-only; polars is used only for mcols.
Masked positions are NA; their underlying data is unspecified.

## Consequences
Callers use `.mask`/`.filled(...)`; tests compare via `.filled(-1)` after
checking the mask.
