# ranges_e0 — E0 control port of IRanges / GRanges

A plain Python port of Bioconductor's `IRanges` and `GRanges` (scope in
`../../TARGET.md`; reference: R 4.6, IRanges 2.47.2, GenomicRanges 1.65.1).
E0 is the **control** experiment: no conformance suite, spec or ADRs. The
only artifacts are the code, a pytest suite with hand-written expectations,
this README and an append-only judgment log, [NOTES.md](NOTES.md).
[PLAN.md](PLAN.md) is the original plan; where it was wrong, NOTES.md says so.

```sh
uv run pytest          # 300 tests (incl. parametrized), ~3 s
uv run --with mypy mypy src --ignore-missing-imports   # type check (clean)
```

```python
from ranges_e0 import IRanges, GRanges, Seqinfo, NO_MATCH

x = IRanges(start=[1, 5, 20], end=[4, 8, 25])
x.reduce().to_tuples()              # [(1, 8), (20, 25)]
x.find_overlaps(IRanges(3, 6))      # Hits, 0-based (query, subject) pairs
g = GRanges(["chr1", "chr1"], IRanges([10, 30], width=5), ["+", "-"],
            seqinfo=Seqinfo(["chr1"], [1000]))
g.promoters(upstream=5, downstream=2).to_tuples()
```

## Design

**Thin classes over numpy kernels.** Every algorithm is a pure function on
int64 `start`/`width` arrays in a private module; the classes validate, call a
kernel and rebuild an object.

| module | content |
|---|---|
| `iranges.py` | `IRanges` (ctor, accessors, subsetting, `concat`, method wrappers), `NO_MATCH` |
| `granges.py` | `GRanges` (seqname codes + strand codes + `Seqinfo`), the grouping seam |
| `seqinfo.py` | `Seqinfo` (+ `merge`) |
| `rle.py`, `hits.py` | minimal `Rle` (coverage result) and `Hits` |
| `_intra.py` | shift / narrow / resize / flank / promoters |
| `_inter.py` | range / reduce / gaps / disjoin |
| `_setops.py` | union / intersect / setdiff and `p*` versions |
| `_overlaps.py` | the overlap predicate and `overlap_pairs` engine |
| `_nearest.py` | precede / follow / nearest / distance |
| `_coverage.py`, `_order.py`, `_validate.py` | coverage sweep, ordering/rank/duplicated, argument coercion |

**GRanges = IRanges kernel per group.** `GRanges._groups()` splits rows by
`(seqlevel, strand)` (or seqlevel only with `ignore_strand`) and `_assemble()`
glues per-group results back in R's output order (seqlevel in `Seqinfo` order,
then `+ - *`). Strand-aware operations add only a small amount of strand logic
on top (flip start/end for `-` in `resize`/`flank`/`promoters`; direction
reversal in `precede`/`follow`; `*` compatibility in overlaps and `p*` ops).

**Overlap engine.** Subjects are sorted by start; each query gets a contiguous
candidate block from `searchsorted` (upper bound from subject starts, lower
bound from the running maximum of subject ends); the exact predicate then
filters candidates, in chunks to bound memory. 200k × 200k random ranges take
about a second. The predicate for each `type`/`maxgap`/`minoverlap` is
documented at the top of `_overlaps.py`; it was inferred from R behaviour
tables because the C code is not available here.

**Storage.** `start` and `width` are int64 (R: int32, so R's overflow errors
are not reproduced); `end = start + width - 1`; zero-width ranges have
`end = start - 1`. `mcols` is a `polars.DataFrame` (zero-column frame = none).
Objects are immutable by convention; R's replacement functions are
`with_names`, `with_mcols`, `with_strand`, `with_seqinfo`.

## API (R name → Python)

Methods use snake_case names and arguments (`min.gapwidth` → `min_gapwidth`).
`c(x, y)` → `concat(x, y)`; `length(x)` → `len(x)`.

| area | Python (IRanges and GRanges unless noted) |
|---|---|
| construct | `IRanges(start, end, width, names, mcols)`; `GRanges(seqnames, ranges, strand, mcols, seqinfo, seqlengths, names)`; `Seqinfo(seqnames, seqlengths, is_circular, genome)` |
| access | `start end width names mcols` · GRanges also `seqnames strand ranges seqinfo seqlevels seqlengths is_circular genome` |
| subset | `x[int]`, `x[slice]`, `x[int array/list]`, `x[bool array]`, `x["name"]`, `x[["a", "b"]]` |
| intra-range | `shift narrow resize(fix=) flank(start=, both=) promoters` (GRanges: `ignore_strand` where R has it) |
| inter-range | `range reduce(min_gapwidth, drop_empty_ranges, with_revmap) gaps(start, end) disjoin(with_revmap)` |
| set ops | `union intersect setdiff punion(fill_gap) pintersect(resolve_empty / GRanges: drop_nohit_ranges, strict_strand) psetdiff` |
| overlaps | `find_overlaps(type, maxgap, minoverlap, select, drop_self, drop_redundant) count_overlaps overlaps_any subset_by_overlaps(invert)` |
| nearest | `nearest(select) precede(select) follow(select) distance` |
| coverage | `coverage(shift, width, weight)` → `Rle` (GRanges: `dict[str, Rle]`) |
| ordering | `order sort rank(ties_method) duplicated unique` (GRanges: `ignore_strand`) |
| seqinfo | `Seqinfo.merge`, `GRanges.trim`, `GRanges.with_seqinfo`; combine merges seqinfo |

## Conventions and deliberate deviations from R

Full reasoning in NOTES.md; the ones a user will hit:

- **Index outputs are 0-based** (`Hits`, `order`, `nearest`/`precede`/`follow`,
  `select=` results, `revmap`). **"No match" is `NO_MATCH == -1`** (also used by
  GRanges `distance` for incompatible pairs). `-1` is a valid Python index, so
  check for it. `rank` stays 1-based because ranks are not indices.
- Negative subscripts follow Python (`x[-1]` is the last element), not R's
  exclusion semantics.
- `Hits` are always sorted by (query, subject). R's within-query subject order
  is an NCList traversal artifact (e.g. 1,3,2,5); compare as sorted pairs.
- `select="arbitrary"` returns the first candidate.
- `rank` defaults to `ties_method="first"` (R's generic has no default and errors).
- `coverage` of a GRanges is a plain dict, not an `RleList`; `Rle` has no arithmetic.
- Only length-1 arguments are recycled (R also recycles other lengths with a warning).
- Out-of-bound ranges warn with `OutOfBoundWarning` (R warns too; text differs).
- Not implemented: GRangesList/IRangesList, Views, `distanceToNearest`,
  `selectNearest`, `overlapsRanges`, integer vectors as query, circular-sequence
  wrap-around, `seqinfo<-` pruning/renaming modes, `RleList`, file I/O.

R quirks that **are** reproduced (and tested): `reduce` keeps isolated empty
ranges by default and merges nearby empties when `min_gapwidth > 1`; `gaps`
default bounds ignore empty ranges and GRanges `gaps` reports every
seqlevel × strand when lengths are known; `disjoin` splits at empty ranges;
`findOverlaps(type="within", maxgap=0)` means "no limit"; `pintersect` on
GRanges adds a `hit` column; GRanges `nearest` resolves ties differently from
IRanges `nearest`; `*`-strand queries consult both directions in
`precede`/`follow`.

## Tests

`tests/` holds ordinary pytest files, one per area (`test_iranges_*.py`,
`test_granges_*.py`, `test_seqinfo.py`, `test_coverage.py`) plus
`test_properties.py` (hypothesis). Two kinds of test:

- **Example tests** with hand-typed expected values. Where a value was checked
  in R, the comment says what R returned (`# R: ...`); that comment is the
  whole provenance trail. There are no golden files, shared cases or runners
  that execute R.
- **Property tests**: set algebra on covered positions, `reduce` idempotence,
  `disjoin`/`gaps` partitions, sort/rank/unique consistency, coverage totals,
  and the overlap engine against an independent explicit-loop reference for
  every `type`/`maxgap`/`minoverlap` (this found one real bug).

Honest limit of the control: expected values encode the implementer's reading
of R. Cases that were not checked ad hoc in R may be wrong in a way the suite
cannot detect.

## Maintaining this

There is deliberately no automation. When upstream behaviour changes:

1. Read the IRanges / GenomicRanges NEWS entry and the help page.
2. Reproduce the new behaviour with a small `Rscript -e` snippet.
3. Edit the kernel (usually one function in a `_*.py` module); GRanges follows
   through the grouping seam unless the change is strand logic in `granges.py`.
4. Update or add an example test with a `# R:` comment, and add a NOTES.md
   entry if it changes a deliberate deviation.
5. Run `uv run pytest`. Nothing flags a divergence that no test covers — that
   is exactly the gap the other experiments try to close.
