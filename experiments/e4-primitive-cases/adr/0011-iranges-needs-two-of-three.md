# ADR-0011: `IRanges` needs two of start/end/width

- Status: accepted. Divergence from R (deliberate).
- R: `IRanges(start=c(1,2,3,7))` is the "integer set" form and returns the
  runs `[1,3]`, `[7,7]`. `IRanges(end=...)` alone is an error.
- Python: fewer than two of start/end/width is always a `ValueError`.
- Why: one rule is easier to state than R's special case; the run-length
  form is `reduce` on width-1 ranges.
- Case: `ir-con-start-only-runs` (`python: {expect_error: true}`).
