# ADR-0006: coordinates beyond int32 raise

- Status: accepted. Divergence from R (deliberate).
- R: `IRanges(start=2147483647L, width=2L)` returns a range whose `end` is
  `NA` plus a warning ("NAs produced by integer overflow").
- Python: `ValueError` at construction (any start/end outside
  [-2^31+1, 2^31-1]).
- Why: NA coordinates silently corrupt every downstream operation.
- Case: `ir-con-int32-overflow` (`python: {expect_error: true}`); the golden
  stays R's output.
