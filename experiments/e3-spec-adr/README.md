# e3-spec-adr — IRanges / GRanges in Python, spec and ADRs first

Experiment e3: a behavioral spec and ADRs were written first (probing R
4.6 / IRanges 2.47.2 / GenomicRanges 1.65.1), then the package and tests
were built against them. What the spec does not describe is
implementation-defined. R is normative unless an ADR says otherwise.

```python
from ranges_e3 import IRanges, GRanges, concat

x = IRanges(start=[1, 5, 10], width=[3, 0, 4])
x.reduce(), x.gaps(start=1, end=20), x.find_overlaps(x.shift(2), type="within")

g = GRanges(["chr1", "chr1"], IRanges(start=[10, 40], end=[20, 50]), strand=["+", "-"])
g.flank(5), g.promoters(100, 10), g.reduce(ignore_strand=True), g.sort()
```

Element indices are 0-based; coordinates are 1-based closed (ADR-0001).
Everything is `uv run pytest` (162 tests).

## Layout

| Path | Role |
|---|---|
| `spec/NN-*.md` | Normative behavior, language-neutral, clause ids `**AREA.n**` |
| `adr/NNNN-*.md` | Non-obvious choices, each linked to clause ids |
| `src/ranges_e3/_kernels.py` | All range math: pure numpy over `(start, width)` int64 arrays |
| `src/ranges_e3/iranges.py`, `granges.py`, `seqinfo.py`, `hits.py`, `rle.py` | Thin classes: validate, call kernels, rebuild objects |
| `tests/test_*.py` | One file per spec area; each test cites clauses with `@pytest.mark.spec("OVL.3")` |
| `tests/test_spec_trace.py` | Fails if a live clause has no citing test, or a test cites an unknown id |
| `NOTES.md` | Timestamped judgment calls (append-only) |
| `PLAN.md` | The original plan (deviations are in NOTES.md) |

Design: deep module, clear seam. `GRanges` only groups rows by
`(seqlevel, strand)` or builds comparable blocks (`_blocks`) and feeds the same
kernels `IRanges` uses; strand flips are passed to kernels as per-range arrays.
`mcols` are polars frames, NA ints are masked arrays, `Hits` and `Rle` are
minimal result types (ADR-0001..0005).

## Spec areas

| File | Prefix | Covers |
|---|---|---|
| 00-conventions | CONV | notation, indices, NA, recycling, names/mcols propagation |
| 10-iranges-core | CORE | construction, accessors, subsetting, concat |
| 20-intra | INTRA | shift, narrow, resize, flank, promoters |
| 30-inter | INTER | range, reduce, gaps, disjoin |
| 40-setops | SET | union/intersect/setdiff and parallel versions |
| 50-coverage | COV | coverage, `Rle` |
| 60-overlaps | OVL | find_overlaps, counts, `Hits` |
| 70-nearest | NEAR | nearest, precede, follow, distance |
| 80-order | ORD | order, sort, rank, duplicated, unique |
| 90-seqinfo | SEQI | Seqinfo, merge |
| 95-granges | GR | GRanges: strand rules, grouping, trim, ordering |

## Maintaining the artifacts

- **Behavior change** (e.g. upstream NEWS): edit the clause in place, keep its
  id; update or add `Ex:` lines; change code; fix the citing tests. If R's new
  behavior is surprising or the port deliberately differs, add an ADR.
  Commit spec and code together or spec first (`spec: amend OVL.4 — why`).
- **Never renumber.** Retire a clause by keeping the id, striking the text and
  adding `Retired: <why>` on the clause line; the trace test then skips it.
- **ADRs are append-only.** Supersede with a new ADR; set the old status to
  `superseded by ADR-NNNN`.
- **New clause ⇒ new test.** Cite the id in `@pytest.mark.spec(...)`.
  `test_spec_trace.py` enforces this; its `PENDING_AREAS` allow-list is for
  stages in progress and is empty now.
- **Spec is not a test suite.** Clauses state rules and a few examples; the
  tests check more cases but never define behavior.
- **Checking against R** is manual: run `Rscript` on the clause's `Ex:` lines.
  No R goldens or cross-language harness are committed. During development a
  gitignored scratch fuzz (`spec/scratch/`) compared ~100k random cases with R;
  it found the corner cases now recorded in INTER.3–6, INTRA.2, GR.11–15 and
  ADR-0013/0014. Re-create it (R emits cases + answers, Python recomputes) when
  touching kernels.

## Known divergences from R

- `Hits` order within a query is subject-index ascending (ADR-0009).
- GRanges `nearest(select="all")` keeps all minimal-distance ties; R drops some (ADR-0011).
- Circular sequences get no wrap-around (ADR-0014).
- Self-`nearest` and `precede`/`follow` on zero-width ranges in R-quirk corners
  are implementation-defined (NEAR.2, GR.15).
- Length mismatches are refused where R warns and recycles (ADR-0006); int64
  instead of int32 (ADR-0005); no `weight` in `coverage`, no `Views`/lists/IO.

## Performance

Correctness over speed. Overlaps use a start-sorted candidate window
(sorted subject starts + maximum subject width), chunked to bound memory; worst
case (one huge subject range) is quadratic. `reduce(drop_empty_ranges=True)`
uses a Python loop (ADR-0007).
