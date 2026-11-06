# ADR-0004: R exclusion subscripts are not ported

- Status: accepted. Divergence from R (deliberate).
- R `x[-1]` drops element 1. Python negative ints follow numpy (`x[-1]` is
  the last element).
- Cases use positive indices only; `0` is never used (R: drop; numpy: last
  after the -1 shift).
