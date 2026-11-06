# E0 control — implementation plan

Plain agentic port of IRanges / GRanges (scope: `../../TARGET.md`) to a
Python package `ranges_e0`. No cross-language process: R is the
reference, consulted through help pages, package source and small ad hoc
`Rscript -e` snippets. Tests are ordinary pytest with expected values
written by hand.

Reference versions on this machine: R 4.6.0, IRanges 2.47.2,
GenomicRanges 1.65.1, S4Vectors 0.51.6 (Seqinfo is now a separate
package).

---

## 1. Rules the implementer must keep (from PROTOCOL / task brief)

- Write only inside `experiments/e0-control/`.
- No BiocPy (`iranges`, `genomicranges`, `biocframe`, `biocutils`): do not
  install, import, read or copy.
- R use: `?topic`, reading R source under the R library, and **small ad
  hoc snippets** to settle what a function does. **Not allowed:** any
  script that runs R and Python side by side, case files shared between
  them, golden/fixture files generated from R output, or bulk-pasting R
  output into tests. Expected values are typed by hand. Checking a value
  in R first and then writing it into a test by hand is fine.
- Commit at every stage below (`git commit`, no squash/amend/rebase).
  Message says what changed and why; end with
  `Co-Authored-By: Claude <noreply@anthropic.com>`.
- Judgment calls go in `NOTES.md` (append-only, timestamped
  `## YYYY-MM-DD HH:MM — topic`).
- Foreground only. `uv run pytest` must pass at every commit.

---

## 2. Architecture

### 2.1 Design principle: thin classes over array kernels

Every algorithm is a **pure function on numpy arrays** (`start`, `width`
or `end`, int64) in a private kernel module. `IRanges` and `GRanges` are
thin wrappers that validate input, call kernels and rebuild objects.
`GRanges` reuses the same kernels by **grouping** rows by
`(seqname, strand)` (or just `seqname` when `ignore_strand=True`, or when
the op is strand-blind) and running the IRanges kernel per group.

This gives one place to get each algorithm right, and GRanges behaviour
comes from "group, then call the IRanges kernel" plus a small amount of
strand logic (flipping direction on `-`).

### 2.2 Data model

**Coordinates.** Keep R's 1-based closed convention: `end = start + width - 1`;
zero-width range has `end = start - 1`. Store `start: int64[n]` and
`width: int64[n]`; `end` is derived. Width must be `>= 0`. int64 rather
than R's int32: we don't reproduce R's integer overflow errors (record in
NOTES).

**`IRanges`** (`iranges.py`)

| field | type | notes |
|---|---|---|
| `_start` | `np.ndarray[int64]` | read-only view exposed as `.start` |
| `_width` | `np.ndarray[int64]` | `.width`; `.end` computed |
| `_names` | `np.ndarray[object] \| None` | str per range; `None` = no names |
| `_mcols` | `polars.DataFrame` | height n, zero columns by default |

Treat objects as immutable: operations return new objects. Provide
`with_names(...)`, `with_mcols(...)` instead of setters (the R `names<-`
/ `mcols<-` equivalents).

Construction: `IRanges(start=None, end=None, width=None, names=None, mcols=None)`.
Exactly two of start/end/width, or all three when consistent (R errors if
inconsistent). Scalars broadcast like R's recycling for length-1 args.
Also `IRanges.empty()`.

**`Rle`** (`rle.py`). Minimal run-length vector: `lengths: int64[]`,
`values: np.ndarray`. Used for `coverage` output (and internally for
seqnames/strand grouping if handy). Support `len`, `__eq__`,
`to_numpy()` (decompress), `runs` iteration, `__repr__`. No arithmetic.
TARGET says Rle as a general class is out of scope; keep it this small.

**`Hits`** (`hits.py`). `query_hits: int64[]`, `subject_hits: int64[]`,
`query_length`, `subject_length`. Sorted by query then subject (R order).
Methods: `len`, iteration over pairs, `to_numpy()` (n×2), `__eq__`.

**`Seqinfo`** (`seqinfo.py`). `seqnames: tuple[str, ...]`,
`seqlengths: np.ndarray[int64]` with a sentinel for NA (see 2.4) ,
`is_circular: np.ndarray[object]` (True/False/None), `genome` per seqlevel
(`str | None`). `Seqinfo.merge(other)` follows R: union of seqlevels
keeping order of first-seen, error on conflicting non-NA lengths
(R message: "sequence chr1 has incompatible seqlengths").

**`GRanges`** (`granges.py`)

| field | type | notes |
|---|---|---|
| `_seqnames` | `np.ndarray[int32]` | codes into `seqinfo.seqnames` |
| `_ranges` | `IRanges` | no names/mcols inside; names/mcols live on GRanges |
| `_strand` | `np.ndarray[int8]` | `+1` = `+`, `-1` = `-`, `0` = `*` |
| `_seqinfo` | `Seqinfo` | |
| `_names`, `_mcols` | as IRanges | |

Public accessors return friendly forms: `seqnames` → `np.ndarray[str]`
(or `Rle` — pick one, record it), `strand` → `np.ndarray[str]` of
`"+"`, `"-"`, `"*"`, `ranges` → `IRanges`. Construction accepts strings
(`"chr1:1-5:+"` parse is optional sugar; nice for tests), seqnames as a
list, `ranges` as IRanges or start/end/width, strand as list or scalar,
`seqinfo` or `seqlengths` dict.

### 2.3 Python API conventions

- Methods on the classes, Python-cased R names: `shift`, `narrow`,
  `resize`, `flank`, `promoters`, `range`, `reduce`, `gaps`, `disjoin`,
  `union`, `intersect`, `setdiff`, `punion`, `pintersect`, `psetdiff`,
  `find_overlaps`, `count_overlaps`, `overlaps_any`,
  `subset_by_overlaps`, `nearest`, `precede`, `follow`, `distance`,
  `coverage`, `order`, `sort`, `rank`, `duplicated`, `unique`, `trim`.
  Plus module-level `concat(*objs)` (R `c()`), also exposed as
  `IRanges.concat` / `GRanges.concat`.
- Argument names are snake_case of R's: `min_gapwidth`,
  `drop_empty_ranges`, `with_revmap`, `ignore_strand`, `maxgap`,
  `minoverlap`, `type`, `select`, `fix`, `use_names`, `both`,
  `upstream`, `downstream`, `start`/`end` for `gaps`. Same defaults as R.
- `len(x)` is R's `length(x)`.
- Subsetting `x[i]`: int, slice, int array, bool array (length n), str or
  list of str (names; first match like R). Negative ints follow
  **Python** semantics (from the end), not R's exclusion — record in NOTES.

### 2.4 Index base and NA — decide once, apply everywhere

- **Index-valued outputs are 0-based** (`Hits`, `order`, `nearest`,
  `precede`, `follow`, `find_overlaps(select="first"/...)`, revmap
  entries). This is the only usable choice for Python subsetting.
- **"No match" in index outputs = `-1`**, exported as
  `ranges_e0.NO_MATCH`. Caveat: `-1` is a valid Python index, so docs and
  tests must make this visible. Alternative (masked arrays / polars
  nullable Series) was judged clumsier; if the implementer finds `-1`
  causing bugs, a switch is a recorded redesign.
- `distance` with no defined value (different seqnames, incompatible
  strand) → also `-1` (distances are otherwise `>= 0`).
- `rank` returns ranks, which are **not** indices: keep R's 1-based ranks.
- Seqlengths NA → `None` at the API surface; internally `-1` sentinel.

### 2.5 Core algorithms

All in kernels; n = number of ranges.

- **intra-range** (`_intra.py`): vectorised arithmetic. `resize`
  fix="center" rounding (R: `IRanges(10,20)` resize 4 center → 13-16);
  `promoters` on width 0 and on `-` strand; `flank` with
  `start`/`both`; `narrow` with negative start/end/width (relative to the
  range, R semantics). GRanges: `-` strand swaps the meaning of
  start/end for `resize`, `flank`, `promoters`; `*` behaves as `+`.
- **ordering** (`_order.py`): `np.lexsort`. IRanges key: start, then
  width. GRanges key: seqlevel code, strand (`+` < `-` < `*`), start,
  width. `order(decreasing=...)` must be stable. `duplicated` on the
  same key tuple; `unique` = keep first. `rank`: R's generic errors if
  `ties.method` is not given explicitly (verified ad hoc); Python
  default `"first"` and support `"first"`, `"min"`, `"max"`,
  `"average"` — record the deviation.
- **reduce** (`_inter.py`): sort by (start, end); sweep, merging when
  `next.start <= cur_end + min_gapwidth`. Zero-width ranges: dropped
  unless `drop_empty_ranges=False` (verified: IRanges(5,4) survives with
  `drop.empty.ranges=FALSE`; check how a zero-width range adjacent to or
  inside a non-empty range is treated). `with_revmap` → polars `List(Int64)`
  mcol `revmap`, 0-based, ascending within each group.
- **range**: min start / max end (per group for GRanges).
- **gaps**: complement of reduce within `[start, end]` (defaults
  min/max of x for IRanges). GRanges: per seqlevel × strand bounded by
  `[1, seqlength]`. Verified quirks: with known seqlengths, gaps returns
  rows for *every* seqlevel and *all three* strands, including empty
  ones; with NA seqlength only `*`-strand gaps up to the last range
  appear. Read `GenomicRanges:::gaps` source before implementing.
- **disjoin**: breakpoints = all starts ∪ all (end+1); segments between
  consecutive breakpoints kept if covered (prefix-sum of +1/-1 events).
  Check zero-width handling in R.
- **set ops** (`_setops.py`): `union = reduce(concat(x, y))`;
  `intersect`/`setdiff` via sweeps on reduced inputs (R defines
  `setdiff(x, y)` via gaps). Results drop names/mcols like R.
  `p*` versions are elementwise with R's rules (`punion` requires overlap
  or `fill_gap=True`; `pintersect` with `resolve_empty`; `psetdiff`
  errors if the result would be two ranges). GRanges `p*`: same seqname
  required, strand compatibility.
- **overlaps** (`_overlaps.py`): one kernel
  `overlap_pairs(q_start, q_end, s_start, s_end, *, type, maxgap, minoverlap) -> (qi, si)`.
  Implementation: sort subject by start; per query, `searchsorted` gives
  the candidate block with `s_start <= q_end + maxgap + 1`, then filter
  with the full predicate in numpy. Keep a brute-force reference
  (`_overlap_pairs_bruteforce`, O(n·m) broadcasting) **in the test suite**
  and property-test the two against each other with hypothesis. Semantics
  to get exactly from `?findOverlaps`: `maxgap` default `-1`,
  `minoverlap` default `0`, how each applies to `start`/`end`/`within`/
  `equal`, and zero-width ranges (verified: zero-width query at 5
  overlaps 4-6). `select`: `all` → `Hits`; `first`/`last` → smallest/
  largest subject index or `NO_MATCH`; `arbitrary` → any hit (test only
  membership). GRanges wrapper: group by seqname, compatible strands
  (`*` matches everything; `+`≠`-`) unless `ignore_strand`, then map back
  to original indices and sort.
- **nearest family** (`_nearest.py`): `precede` = subject range that x
  precedes (closest one starting after x ends); `follow` mirror;
  `nearest` = overlapping first (any overlap counts as distance 0), else
  closest of precede/follow, `select="arbitrary"|"all"` (ties →
  all; verified: IRanges(5,6) vs 1-3 and 8-9 returns both). `distance`
  is parallel: `0` for adjacent and overlapping (verified
  distance(1-3, 4-5) = 0), else gap width... check the exact definition
  (it changed historically). GRanges: strand-aware — on `-` strand
  "precede" looks to lower coordinates (verified: chr1:1-2:- precede
  chr1:5-6:- → NA); `*` rules need checking in `?precede`.
- **coverage** (`_coverage.py`): event sweep (+1 at start, -1 at end+1),
  cumulative sum, run-length encode → `Rle`; `width` argument sets
  length (default max end); optional `shift`/`weight` only if cheap.
  GRanges → `dict[str, Rle]` over all seqlevels, length = seqlength or
  max end.
- **trim**: clip to `[1, seqlength]` for non-circular seqlevels with
  known length. GRanges construction / shift producing out-of-bound
  ranges emits `warnings.warn` like R (R warns, doesn't error).
- **concat**: IRanges mcols are aligned by column name, missing columns
  filled with null (verified: `c()` fills `<NA>`). GRanges also merges
  seqinfo (error on conflicting lengths, verified).

### 2.6 Module layout

```
experiments/e0-control/
  pyproject.toml           uv project, deps numpy + polars; dev pytest + hypothesis
  PLAN.md                  this file
  NOTES.md                 append-only judgment log
  README.md                design + maintenance (filled in at the end)
  src/ranges_e0/
    __init__.py            public re-exports, NO_MATCH, concat
    rle.py                 Rle
    hits.py                Hits
    iranges.py             IRanges (data model, ctor, accessors, subset, concat, methods)
    seqinfo.py             Seqinfo
    granges.py             GRanges (+ grouping helper)
    _validate.py           arg coercion: int arrays, recycling, names
    _intra.py              shift/narrow/resize/flank/promoters kernels
    _order.py              order/rank/duplicated kernels
    _inter.py              range/reduce/gaps/disjoin kernels
    _setops.py             union/intersect/setdiff + p* kernels
    _overlaps.py           overlap_pairs kernel + select handling
    _nearest.py            precede/follow/nearest/distance kernels
    _coverage.py           coverage kernel
  tests/
    test_iranges_core.py   ctor/accessors/subset/concat/mcols
    test_iranges_intra.py
    test_iranges_order.py
    test_iranges_inter.py
    test_iranges_setops.py
    test_iranges_overlaps.py
    test_iranges_nearest.py
    test_coverage.py
    test_seqinfo.py
    test_granges_core.py
    test_granges_intra.py
    test_granges_order.py
    test_granges_inter.py
    test_granges_overlaps.py
    test_properties.py     hypothesis invariants + kernel-vs-bruteforce
```

Keep files merged if one turns out tiny; don't create modules for the
sake of the layout.

---

## 3. Artifacts and how they are maintained

E0 deliberately has no extra process artifacts. What exists:

1. **pytest suite** — the only executable statement of expected
   behaviour. Two kinds of tests:
   - *Example tests*: small inputs, hand-written expected values. Where a
     value was checked in R, add a short comment such as
     `# R: reduce(IRanges(c(1,5,20), c(2,6,20)), min.gapwidth=3) -> 1-6, 20-20`.
     That comment is the whole provenance trail.
   - *Property tests* (hypothesis): invariants that need no R —
     `reduce` idempotent and sorted/non-overlapping; `union`/`intersect`/
     `setdiff` agree with set algebra on the covered integer positions;
     `coverage` total = sum of widths; `disjoin` pieces don't overlap and
     cover the same positions as `reduce`; `gaps` ∪ `reduce` = full span;
     `sort(x) == x[order(x)]`; `count_overlaps == bincount(find_overlaps)`;
     overlap kernel == brute force for every `type`/`maxgap`/`minoverlap`.
2. **NOTES.md** — append-only, timestamped log of judgment calls and
   deliberate deviations from R (index base, NA sentinel, rank default,
   negative indexing, int64, anything unimplemented).
3. **README.md** — design summary (data model, kernels/grouping seam,
   API conventions, deviations list pointing at NOTES) and a short
   "maintaining this" section: when upstream R behaviour changes, a
   developer reads NEWS, checks it ad hoc in R, edits the kernel, updates
   the example tests. No automated detection — that is the point of the
   control.
4. **Docstrings** on public methods: one line plus R-name mapping and
   any deviation.

---

## 4. Staged work breakdown (one or more commits per stage)

Each stage: implement, write tests, `uv run pytest`, commit. Use ad hoc R
to check anything you're unsure about **before** writing the expected
value.

| # | Stage | Content | Commit message gist |
|---|---|---|---|
| 0 | Scaffold | *done by planner*: uv project, empty package, smoke test, PLAN.md, NOTES.md | — |
| 1 | IRanges core | `_validate`, `IRanges` ctor (any 2 of 3, consistency check, recycling), accessors, names, `__getitem__` all index kinds, `concat` with mcol alignment, `with_names`/`with_mcols`, `__eq__`, `__repr__` (R-like table) ; `Rle`, `Hits` skeletons | "IRanges data model: int64 start/width + polars mcols" |
| 2 | Intra-range | `shift`, `narrow`, `resize` (start/end/center), `flank` (start, both), `promoters` | "IRanges intra-range ops" |
| 3 | Ordering | `order` (decreasing), `sort`, `rank` (ties), `duplicated`, `unique` | |
| 4 | Inter-range | `range`, `reduce` (min_gapwidth, drop_empty_ranges, with_revmap), `gaps` (start/end), `disjoin` | |
| 5 | Set ops | `union`, `intersect`, `setdiff`, `punion`, `pintersect`, `psetdiff` | |
| 6 | Overlaps | `overlap_pairs` kernel, all `type`s, `maxgap`, `minoverlap`, `select`; `find_overlaps`, `count_overlaps`, `overlaps_any`, `subset_by_overlaps`; brute-force property test. Self-overlap form `x.find_overlaps()` (subject omitted) only if R's semantics are clear | "overlap engine: sorted-start + searchsorted ..." |
| 7 | Nearest | `precede`, `follow`, `nearest` (select arbitrary/all), `distance` | |
| 8 | Coverage | `coverage` → `Rle` | |
| 9 | Seqinfo + GRanges core | `Seqinfo` (+ merge), `GRanges` ctor, accessors (`seqnames`, `strand`, `ranges`, `seqinfo`, `seqlevels`, `seqlengths`, `is_circular`, `genome`), subset, `concat` with seqinfo merge, `trim`, out-of-bound warning | |
| 10 | GRanges intra | strand-aware `resize`/`flank`/`promoters`, plain `shift`/`narrow`, `ignore_strand` where R has it | |
| 11 | GRanges ordering | key seqlevel, strand, start, width; `order/sort/rank/duplicated/unique` | |
| 12 | GRanges inter + set ops | grouping helper; `range`, `reduce`, `gaps` (seqlengths/NA quirks), `disjoin`, set ops + `p*`, all with `ignore_strand` | |
| 13 | GRanges overlaps, nearest, coverage | strand compatibility, `ignore_strand`, cross-seqname NA in `distance`; coverage per seqlevel | |
| 14 | Hardening + docs | fill property tests, sweep TARGET.md scope list item by item (tick each to a test), README design + maintenance section, final NOTES entry | "README: design + maintenance; scope audit" |

Ordering rationale: overlap and nearest depend on ordering/sorting
helpers; GRanges depends on every IRanges kernel; set ops depend on
reduce/gaps. If the ~4 h ceiling approaches, prefer **breadth** (every
op public and tested on the basic case) over edge-case depth, and log
what was skipped in NOTES.

---

## 5. Main risks

1. **Zero-width ranges.** They behave specially in reduce, disjoin,
   overlaps (zero-width at 5 overlaps 4-6), promoters and coverage. Test
   each op on a zero-width input and check R ad hoc.
2. **`findOverlaps` parameter interplay.** `maxgap`/`minoverlap` mean
   different things per `type`; `minoverlap > 0` with `maxgap >= 0` is an
   R error in some cases. Read `?findOverlaps` "Details" fully before
   coding the predicate.
3. **Strand semantics in GRanges.** Which ops flip on `-`, how `*` is
   treated in precede/follow/nearest, and how `ignore_strand` changes
   result strand (R sets `*`). Easy to get plausible-but-wrong.
4. **`gaps` on GRanges** produces rows for empty seqlevel/strand
   combinations depending on whether seqlengths are known.
5. **Index base / NA sentinel.** `-1` silently indexes the last element.
   Mitigation: one exported constant, tests that assert `NO_MATCH`
   explicitly, docstrings.
6. **Hand-written expectations can encode the implementer's
   misunderstanding** — the control has no systematic R check. The
   mitigation within E0 rules is ad hoc R checks for every non-obvious
   case and the `# R:` comment habit. This is the risk the experiment is
   designed to measure; don't try to "fix" it by building a harness.
7. **R-ordering of outputs** (Hits sorted by query then subject; reduce/
   gaps output ordered by seqlevel then strand `+ - *`; revmap order).
   Tests should assert order, not just set equality.
8. **Seqinfo merge** order of seqlevels and conflict rules.
9. **Scope creep**: Rle arithmetic, GRangesList, `findOverlaps` on
   self with `drop.self`/`drop.redundant`, `weight` in coverage. Skip
   unless trivial.
10. **Time ceiling** (~4 h): stages 1–8 must not eat the whole budget;
    timebox each IRanges stage.
