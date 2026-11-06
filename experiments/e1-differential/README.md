# ranges_e1 — IRanges / GRanges in Python (experiment e1-differential)

Independent numpy/polars port of Bioconductor's IRanges 2.47 / GenomicRanges 1.65.
**R is normative.** A differential harness runs every case in R and in Python and
compares the results. Divergences are fixed in Python or recorded in a short ADR.

```
uv run pytest                                   # 3.8k tests, ~10 s, no R needed (R results are cached)
uv run pytest tests/differential -k reduce      # a subset
uv run pytest --refresh-r                       # rerun everything in R, rewrite the cache
uv run python tests/differential/show.py <id>   # R vs Python side by side for one case
```

## Design

| Module | Role |
|---|---|
| `iranges.py` | `IRanges`: immutable `(start, width)` int64 arrays + names + mcols |
| `granges.py` | `GRanges`: seqlevel codes + `IRanges` coordinates + strand codes + `Seqinfo` |
| `seqinfo.py` | `Seqinfo` (frozen) with R's merge rules |
| `hits.py`, `rle.py` | `Hits` (0-based index pairs), minimal `Rle` for coverage |
| `_kernels/` | pure numpy, no classes: `intra inter setops overlaps nearest coverage order` |
| `_meta.py` | names / mcols / subscript helpers |

**Seam.** Kernels take flat arrays plus a group key. `IRanges` passes one group;
`GRanges` passes `seqlevel*3 + strand` (or just the seqlevel with `ignore_strand`) and
mirrors strand-sensitive arguments for `-`. Interval logic lives once, in `_kernels/`.

**Conventions** (ADR-0001/0002)
- Coordinates are 1-based closed like R. *Positional* results (Hits, `order`, `nearest`,
  `revmap`, `select=`) are **0-based**; NA is a masked array entry, never `-1`.
- R generics → snake_case methods (`findOverlaps` → `find_overlaps`, `min.gapwidth` →
  `min_gapwidth`). `x.reduce()`, `x.intersect(y)`, `IRanges.concat([...])`.
- Objects are immutable; `names<-`, `mcols<-`, `seqlengths<-`… are `set_names`, `set_mcols`,
  `set_seqlengths`… returning copies. `mcols` is a `polars.DataFrame` or `None`.
- R warnings → `UserWarning`; R errors → `ValueError`/`IndexError`.
- Vector results are plain arrays (no R-style names, ADR-0003); `seqlengths`, `is_circular`,
  `genome` are dicts keyed by seqlevel; `coverage` is an `Rle` (IRanges) or `dict[str, Rle]` (GRanges).
- O(n·m) overlap search (chunked broadcasting). Correctness first; not tuned.

Public surface: everything in TARGET.md — see the method tables in `iranges.py` / `granges.py`
(every operation has ≥1 differential case).

## Differential harness

Case → `{inputs, call}` in **R vocabulary** (YAML) → R evaluates it (`r/run_cases.R`),
Python evaluates it through `adapter.py` → both emit **canonical JSON** (R's 1-based indices) →
`codec.first_diff` compares. Full syntax and how to add a case: `tests/differential/README.md`.

- `cases/*.yaml`: handwritten, one file per op group; `cases/grid_*.yaml` are **generated** by
  `tests/differential/gen/*.py` (seeded; exhaustive small pools or random pools). Edit the generator.
- `r_cache/*.json`: committed R results, keyed per case by sha256 of its input/call. The whole file is
  re-run when R/IRanges/GenomicRanges versions or the R harness scripts change; stale + no R → skip.
- Python-internal errors (`TypeError`, `AttributeError`…) fail a case even when "both sides error".
- `divergence: {adr, python}` pins an ADR-sanctioned Python value so it stays a regression test.

### Maintaining the artifacts
| Artifact | Rule |
|---|---|
| `cases/*.yaml` | append; edit only to fix a mistranslated case (say so in the commit) |
| `r_cache/` | harness-generated, never hand-edited; commit with the case that produced it |
| `adrs/` | only for accepted divergences / real design choices, linked from cases via `divergence.adr` |
| `NOTES.md` | append-only, timestamped judgment calls and R quirks found |
| Maintenance replay | upgrade R packages → `uv run pytest` reruns R (version key changed) → `git diff r_cache/` is the list of upstream behaviour changes → fix Python or add a case/ADR |

### Divergences from R (ADRs)
| ADR | Python differs from R |
|---|---|
| 0003 | `countOverlaps` returns a plain array (R carries query names) |
| 0004 | `rank()` defaults to `ties_method="average"` (R's default errors) |
| 0005 | int32 overflow raises or is computed exactly (R wraps to NA / invalid objects) |
| 0006 | `mcols()` works when names contain NA (R errors) |

### R quirks reproduced (not ADRs; see NOTES.md for evidence)
`reduce(min.gapwidth=0)` drops repeated empty ranges; `reduce(drop.empty.ranges=TRUE)` lets empties
bridge/extend clusters; `disjoin` splits at empty ranges but never reports them; hit order and
`select="arbitrary"` depend on which side R builds its NCList on (per seqlevel for GRanges);
GRanges `nearest` ties go to the lower index, IRanges ties to `follow`; GRanges precede/follow
drop subjects with coordinate 0 on non-first seqlevels (sentinel collision).

## Not implemented
`IRanges(start=)`-style coercions beyond integer/logical runs, multi-argument `order(x, y)`,
`rank(ties_method="random")`, `seqlevels<-` (renaming/pruning), `GRangesList`/`IRangesList`/Views,
file I/O, `mcols` row names, `circle.length` arguments, comparison operators (`==`, `<`, `match`).

## Status
3720 differential cases and 220 unit tests: **3888 passed, 52 xfailed** after red-team testing.
The 119 new cases expose 47 bug cases and 5 ADR judgments; differences remain unresolved.
[Red-team results](REDTEAM.md) link every case's R/Python comparison. `src/` was unchanged.
The implementation-stage `mypy src` check was clean (`uvx mypy src --ignore-missing-imports`).

### Where the harness has been pushed (and where to push next)
Generators in `tests/differential/gen/` cover: exhaustive small pools (parallel ops, overlaps, nearest),
seeded random pools (IRanges/GRanges, large nested inputs), metadata propagation, mcols type mixes,
malformed arguments to every op, constructors, and random seqinfo pairs (warnings compared).
The red team added circular-origin/periodic overlaps, empty/zero-width interactions, ordering ties,
metadata recycling, missing weights, and validation adversaries. See `REDTEAM.md`.
Untested by design: values near ±2^31 (ADR-0005), `random` ties, `seqlevels<-`, performance at 1e5+ ranges.
