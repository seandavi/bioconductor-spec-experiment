# ADR-0009: `count_overlaps` returns unnamed counts

- Status: accepted. Divergence from R (deliberate).
- R names the integer result by the query's names.
- Python returns a plain int array. Names are a property of the query; the
  caller has `x.names`.
- Case `ir-count-names`: `compare: {fields: {names: ignore}}`.
