# ranges-e0b

IRanges / GRanges in Python, written from knowledge of Bioconductor semantics only
(experiment **e0b-no-r**: no R access, no BiocPy, no conformance harness or formal spec).
Scope: `TARGET.md`. Original plan: [PLAN.md](PLAN.md). Judgment calls: [NOTES.md](NOTES.md).

```sh
uv sync
uv run pytest        # ~190 tests, ~2 s
```

## Honest status

- All TARGET operations are public and tested.
- Expected values in tests are hand-derived. Anything about R I was unsure of is tagged
  `UNVERIFIED` in tests and NOTES.md; nothing was checked against R.

## Design

- **Deep kernels, thin classes.** Algorithms work on numpy arrays `(start, end, group)`.
  `IRanges` uses one group; `GRanges` uses `group = seqlevel_code*3 + strand_code`. One kernel serves both.
- **Immutable values.** Every op returns a new object. Arrays are read-only; mcols frames are copied in
  and out (mutating a frame you passed or got back never changes the object). No setters
  (`with_names`, `with_mcols`, `with_seqinfo`).
- **Vectorised.** Loops only over seqlevels / strands.
- `GRanges` stores seqnames as codes into `seqinfo`, so seqlevels order (not alphabet) drives sorting.
  Strand codes `+`=0 `-`=1 `*`=2 (R factor order).
- mcols are `polars.DataFrame`.

```
src/ranges_e0b/
  _sew.py       start/end/width solving, index resolution, mcols helpers, NO_HIT
  iranges.py    IRanges            granges.py   GRanges        seqinfo.py  Seqinfo
  intra.py      shift/narrow/resize/flank/promoters kernels
  inter.py      reduce/range/gaps/disjoin kernels (grouped)
  setops.py     union/intersect/setdiff sweep + parallel p* kernels
  overlaps.py   Hits + overlap kernel (sorted subject, searchsorted window, exact filter)
  nearest.py    precede/follow/nearest/distance (seq/strand-aware driver)
  coverage.py   rle.py   ordering.py
```

## API conventions (deliberate differences from R)

| Topic | Python | R |
|---|---|---|
| Coordinates | 1-based, closed | same |
| Positional indices (subsetting, `Hits`, revmap, `order`, nearest results) | 0-based | 1-based |
| No hit (`select=`, nearest/precede/follow) | `NO_HIT` = -1 | `NA` |
| Negative subscript | from the end | exclusion |
| Names | snake_case: `min_gapwidth`, `drop_empty_ranges`, `with_revmap`, `ignore_strand`, `maxgap`, `minoverlap` | dotted |
| Ops | methods: `x.reduce()`, `x.find_overlaps(y)` | generics |
| Integers | int64, no overflow check | int32 |
| `reduce(with_revmap=True)` | returns `(ranges, list[ndarray])` | `revmap` mcols column |
| `rank` | 1-based values, float64 for "average" | same |
| Equality | `x.equals(y)` (no element-wise `==`) | `==` |
| Concat | `x.concat(y, z)`, `IRanges.combine([..])` | `c()` |

Results of `reduce/range/gaps/disjoin` and set ops drop names and mcols. Intra-range ops,
subsetting, concat, sort keep them. Parallel ops keep x's.

## R name -> method

| R | Python (IRanges and GRanges) |
|---|---|
| `start end width names length mcols` | `.start .end .width .names len() .mcols` |
| `seqnames strand ranges seqinfo seqlevels seqlengths isCircular genome` | `.seqnames .strand .ranges .seqinfo .seqlevels .seqlengths .is_circular .genome` |
| `shift narrow resize flank promoters` | `.shift .narrow .resize .flank .promoters` |
| `range reduce gaps disjoin` | `.range .reduce .gaps .disjoin` |
| `union intersect setdiff punion pintersect psetdiff` | `.union .intersect .setdiff .punion .pintersect .psetdiff` |
| `findOverlaps countOverlaps overlapsAny subsetByOverlaps` | `.find_overlaps .count_overlaps .overlaps_any .subset_by_overlaps` |
| `nearest precede follow distance` | `.nearest .precede .follow .distance` |
| `coverage` | `.coverage` (IRanges -> `Rle`; GRanges -> `dict[str, Rle]`) |
| `order sort rank duplicated unique` | `.order .sort .rank .duplicated .unique` |
| `trim` | `GRanges.trim` |
| `merge(seqinfo)` | `Seqinfo.merge`, applied automatically on combine and on binary ops |

## Semantics worth knowing

- Empty ranges: `width=0`, `end = start-1`. `reduce` absorbs them into neighbours; `disjoin` never outputs them.
- `findOverlaps` predicate: `gap = max(starts) - min(ends) - 1`; `type="any"` hits iff `gap <= maxgap`
  (default -1 = share a position, adjacency needs `maxgap=0`). Details and empty-range rule in NOTES.
- GRanges `gaps` yields all of `+ - *` per seqlevel (full-length gap for an empty strand).
- Overlaps and nearest on GRanges require the same seqlevel and compatible strand (equal or `*`).
  `type="start"/"end"` flips to the other end for `-`. Circular wraparound is not implemented.
- Seqinfo merge: `None` yields to a value, conflicting values raise, no shared names warns.
- Out-of-bound ranges on non-circular sequences with known length warn (not error); `trim()` fixes them.

## Tests and how to maintain them

- One file per module in `tests/`. Example tests carry an `# R:` comment with the equivalent R call,
  so they can be replayed in R later. Expected values were hand-derived; doubtful ones say `UNVERIFIED`.
- Hypothesis property tests compare against brute force: position sets for reduce/gaps/disjoin/set ops,
  a per-pair oracle for overlaps/nearest/precede/follow, a dense array for coverage, and
  "GRanges op == IRanges op per (seqlevel, strand)".
- When R behaviour is learned, fix the code, replace the `UNVERIFIED` tag in the test, and append to NOTES.md.
- `PLAN.md` is the historical plan; deviations are in NOTES.md, not in PLAN.md.

## Known limits

- Overlap kernel: candidate window per query assumes modest subject widths; one huge subject range
  against many queries degrades to O(n*m). An NCList is the upgrade path.
- No circular wraparound, GRangesList/IRangesList, Views, file I/O.
- Scalar-or-vector arguments are broadcast with length 1 or n only (no R-style full recycling).
