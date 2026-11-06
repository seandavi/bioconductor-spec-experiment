# ADR-0001: Python indices are 0-based

- Status: accepted. Divergence from R (deliberate).
- Scope: every position/index in the Python API: subscripts, `Hits`
  query/subject, `order`, `nearest`/`precede`/`follow`, `select=` results,
  `revmap`.
- Coordinates (`start`/`end`) stay 1-based closed, as in R.
- Goldens stay R-native (1-based). The Python runner adds 1 when projecting
  index fields (`ops.py` `index` flag; `revmap` mcols; `Hits`) and subtracts 1
  from `int` subscript inputs of `[`.
- Why: Python subscripts are 0-based; a 1-based index API would break
  `x[order]`-style idioms. Projection keeps the cases R-faithful.
- `rank` returns rank values (1 = smallest), not positions: unchanged.
