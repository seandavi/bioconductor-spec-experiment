# ADR-0001: 0-based positional results, 1-based coordinates

Status: accepted (plan stage)

## Context
R returns 1-based indices (Hits, `order`, `select="first"`, `nearest`, revmap)
and uses NA for "no match". Python indexes from 0; numpy ints have no NA.

## Options
1. Mirror R: 1-based indices. Every `x[hits]` in Python is off by one.
2. 0-based indices, `-1` for no match. `x[-1]` silently returns the last element.
3. 0-based indices, `numpy.ma.MaskedArray` for NA-able results.

## Decision
Option 3. Genomic coordinates stay 1-based closed (they are data, not
indices). The differential adapter adds 1 to index-valued Python results
before comparing, so the canonical form stays R's.
