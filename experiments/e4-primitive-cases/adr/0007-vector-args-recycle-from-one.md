# ADR-0007: vector arguments recycle only from length 1

- Status: accepted. Divergence from R (deliberate).
- R recycles shorter vector arguments cyclically (`shift(x, c(1,2))` on
  4 ranges works; non-multiples only warn).
- Python: length-1 recycles; any other length != n raises `ValueError`.
- Why: silent cyclic recycling hides length bugs.
- Case: `ir-shift-cyclic` (`python: {expect_error: true}`); golden keeps R's output.
