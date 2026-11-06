# ADR-0008: `seqnames()` / `strand()` project as plain vectors

- Status: accepted. Divergence from R (deliberate), scope decision.
- R returns `Rle` for `seqnames(x)` and `strand(x)`. TARGET.md puts `Rle` as
  a general class out of scope.
- Python returns plain `str` arrays. The R runner expands these two results
  with `as.character()` before projecting (`run_cases.R`, `run_call`), so
  the goldens hold chr vectors and no case needs a compare override.
- `coverage` is the only place `Rle` appears; it keeps run values/lengths.
- Similarly named vectors (`seqlengths`, `isCircular`, `genome`) are Python
  `dict`s (seqname -> value); the projector turns them into values + names.
