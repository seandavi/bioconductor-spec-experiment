# ADR-0005: int32 overflow never wraps

Status: accepted (divergence from R)

## Context
R does coordinate arithmetic in 32-bit integers. Near ±2^31 it silently wraps:
- `shift(<GRanges>, 2147483647)` returns a GRanges whose `start()` is NA (case `gr-shift-overflow-error`);
- `resize(x, 2147483647)`, `flank(x, 2147483647)`, `promoters(x, ., 2147483647)` return objects whose
  `end()` exceeds int32; `findOverlaps(maxgap = 2147483647)` can return no hits; `coverage(shift=)` loses runs.
(IRanges `shift`/constructors do raise.) Probed with `generated bad-args cases` before the value was removed from
`gen/grid_bad_args.py`.

## Options
1. Mirror R: emulate 32-bit wrap-around.
2. Python stores int64: raise `ValueError` when a *result* coordinate leaves int32, and compute
   parameters such as `maxgap` exactly.

## Decision
Option 2. Wrapped results are invalid objects (NA coordinates) that violate the class invariant,
and emulating them would need int32 arithmetic in every kernel. Only the GRanges `shift` case is
pinned with a `divergence:` entry; the other wrap cases are not carried as cases because they pin R's
corrupt output rather than a behaviour anyone should rely on.
