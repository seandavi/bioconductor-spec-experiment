# 00 Conventions (CONV)

Area: notation and rules shared by every other file. Verified against
R 4.6.0, IRanges 2.47.2, GenomicRanges 1.65.1, S4Vectors 0.51.6. NEAR.1 (`distance`) was re-verified against IRanges 2.47.3 / GenomicRanges 1.65.2 behavior (see NOTES 2026-11-06); everything else is as probed on the versions above.

Anything the spec does not describe (error text, repr, performance,
unlisted argument combinations) is implementation-defined.

## Notation

- Range `[s,e]`: closed, 1-based. `w = e - s + 1`. Zero-width: `[s,s-1]`.
- Vector `(1,2,NA)`; set of ranges `{[1,3],[5,4]}`; `x[i]` = i-th element, 0-based.
- `gap(a,b) = max(a.s - b.e, b.s - a.e) - 1`. Negative = overlap, 0 = adjacent.
- "refuse" = the operation raises an error (`ValueError`, or `IndexError` for subscripts).
- Argument names are the Python spelling (`min_gapwidth`).

## Clauses

**CONV.1** Coordinates are 1-based and closed: `end = start + width - 1`. Starts and ends may be zero or negative.

**CONV.2** Width is never negative. Zero-width ranges are legal values; each operation states how it treats them.

**CONV.3** Element indices (subscripts, `Hits`, `nearest`/`precede`/`follow` results, `revmap`, `order`) are 0-based. Rank values are 1-based counts. ADR-0001.

**CONV.4** An integer result that can be missing (`select != "all"` results, `nearest`, `precede`, `follow`, `distance` on `GRanges`) is a masked int64 array; masked = NA. ADR-0002.

**CONV.5** Objects are immutable by convention. Every operation returns a new object; inputs are never changed.

**CONV.6** A per-range argument (`shift`, `width`, `start`, ...) is a scalar (applies to all) or a sequence of length `n` (or length 1, treated as scalar). Any other length is refused. ADR-0006.

**CONV.7** Coordinates are 64-bit integers; no overflow emulation. ADR-0005.

**CONV.8** `names` and `mcols` propagation:

| Operation class | names | mcols |
|---|---|---|
| subset, concat, sort, unique | kept (per CORE.10/CORE.12) | kept |
| shift, narrow, resize, flank, promoters, trim | kept | kept |
| range, reduce, gaps, disjoin, union, intersect, setdiff | dropped | dropped (`reduce`/`disjoin` may add `revmap`) |
| punion, pintersect, psetdiff | from `x` | dropped |
| coverage | n/a | n/a |

**CONV.9** Naming: R `findOverlaps` → `find_overlaps`, `min.gapwidth` → `min_gapwidth`, `ignore.strand` → `ignore_strand`, etc. Two-object operations are methods on the first operand; `concat(*objs)` is a module function. ADR-0004.
