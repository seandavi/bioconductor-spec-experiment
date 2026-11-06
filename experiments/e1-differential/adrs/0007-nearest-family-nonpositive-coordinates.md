# ADR-0007: precede / follow / nearest on GRanges with coordinates <= 0

Status: accepted (divergence from R)

## Context
R's GRanges `precede`/`follow`/`nearest` shift every seqlevel by a multiple of `maxend` and use
`k * maxend` values as sentinels between levels. A range with start/end at or below 0 can collide
with, or cross, those sentinels. Observed in R (cases `rt-edge-gr-negative-*`, found by the red team):
- `nearest(select="all")` returns a hit between *different* seqnames (query 1 on chr2 to subject 3 on chr1);
- `precede(select="all")` returns no hits where a geometric neighbour exists;
- `follow(select="all")` fails with "logical subscript contains NAs";
- an empty range at position 1 (end 0) on the first seqlevel can return a hit on another sequence.

## Options
1. Reproduce the sentinel arithmetic, including cross-sequence hits and the error.
2. Return the geometric answer (same seqname, strand rules unchanged).

## Decision
Option 2. A hit between different sequences is meaningless, and an error on valid input is a defect,
so none of it is worth porting. One narrower quirk *is* kept because ordinary data triggers it: a
subject coordinate of exactly 0 on a seqlevel other than the first is dropped (empty ranges at position 1);
`granges.py::_geometry_hits` documents it. The pinned cases (`divergence:` entries) assert Python's geometric
answers so a future change cannot silently alter them. Generators keep coordinates >= 2 for these operations.
