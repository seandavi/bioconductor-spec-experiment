# ADR-0005: `select="arbitrary"` is compared with `oneof`

- Status: accepted, no divergence in behaviour; looser comparison.
- R returns "an arbitrary hit" per query (observed: the last subject).
  Python returns the lowest. Both are valid.
- Case files use `compare: {fields: {values: oneof}}` plus
  `oneof_from: {kwargs: {select: all}}`. R re-runs the call with
  `select="all"` and stores, per query, the allowed subject indices (or `NA`
  when there is no hit) as `candidates` in the golden.
- Python passes iff each element is in its candidate set.
- Same mechanism applies to `nearest(select="arbitrary")` ties.
