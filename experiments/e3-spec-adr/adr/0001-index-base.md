# ADR-0001 Element indices are 0-based; coordinates stay 1-based

Status: accepted · Clauses: CONV.1, CONV.3, CORE.9, ORD.4

## Context
R uses 1-based coordinates and 1-based element indices. A Python API with
1-based element indices would break `x[i]`, `range`, numpy indexing.

## Options
1. Everything 1-based (R-faithful; fights Python).
2. Everything 0-based, half-open (BiocPy-style; changes every formula).
3. Coordinates 1-based closed (genomic convention, R-faithful); element indices 0-based.

## Decision
Option 3. Subscripts, `Hits`, `nearest`/`precede`/`follow`, `order`, `revmap`
are 0-based; negative subscripts count from the end (Python), not "exclude"
(R). `x[i]` with an int returns a length-1 object. Rank values are counts,
so stay 1-based.

## Consequences
Spec examples need an index-base note (CONV.3). Comparing with R output
means subtracting 1 from element indices only, never from coordinates.
