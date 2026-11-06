# PLAN — e1-differential

Port IRanges/GRanges to Python (`ranges_e1`). R is normative. A
differential harness runs each case in R and Python and compares.
A red-team agent adds cases later. Divergences: fix Python, or write an ADR.

Reference versions: R 4.6.0, IRanges 2.47.2, GenomicRanges 1.65.1,
S4Vectors 0.51.6. `jsonlite` is installed.

## 1. Layout

```
experiments/e1-differential/
  PLAN.md  README.md  NOTES.md  pyproject.toml  uv.lock
  adrs/NNNN-slug.md            # only real decisions
  src/ranges_e1/
    __init__.py                # exports IRanges, GRanges, Seqinfo, Hits, Rle
    _kernels/                  # pure numpy, no classes
      intra.py                 # shift/narrow/resize/flank/promoters on (start, width[, strand])
      inter.py                 # range/reduce/gaps/disjoin per group
      setops.py                # union/intersect/setdiff + p* versions
      overlaps.py              # candidate search + type/maxgap/minoverlap filter
      nearest.py               # nearest/precede/follow/distance
      coverage.py
      order.py                 # lexsort keys, rank, duplicated
    rle.py                     # minimal Rle (values, lengths)
    hits.py                    # Hits
    iranges.py                 # IRanges class: validation, accessors, subset, concat, mcols
    seqinfo.py                 # Seqinfo + merge
    granges.py                 # GRanges class: grouping + strand transforms over kernels
  tests/
    unit/                      # Python-only behaviour (immutability, dtypes, adapter)
    differential/
      README.md                # how to add a case (red-team entry point)
      conftest.py              # --refresh-r option, case collection
      test_differential.py     # one pytest item per case
      harness.py               # load cases, run R (cached), run Python, compare
      adapter.py               # R call vocabulary -> Python API; Python encoder
      codec.py                 # input literal decoding; canonical-value compare
      show.py                  # `uv run python tests/differential/show.py <id>`: R vs Py side by side
      r/run_cases.R            # batch evaluator: JSON cases in, JSON results out
      r/encode.R               # R value -> canonical JSON
      cases/*.yaml             # one file per op group
      r_cache/*.json           # committed R results, keyed by hash
```

**Seam:** kernels take flat arrays plus an optional integer group key.
IRanges calls them with one group. GRanges passes group = seqname code
(× strand code unless `ignore_strand`), and flips strand-sensitive
arguments for `-`. All interval logic lives once, in `_kernels/`.

## 2. Data model

| Type | Fields | Notes |
|---|---|---|
| `IRanges` | `_start: int64[n]`, `_width: int64[n]`, `_names: tuple[str\|None,...] \| None`, `_mcols: pl.DataFrame \| None` | Immutable. Coordinates 1-based closed, `end = start + width - 1`. Validate `width >= 0` and int32 range (R limit). `mcols -> pl.DataFrame \| None`; `None` = zero columns (a zero-column polars frame has height 0, so it cannot carry length). |
| `Seqinfo` | `seqnames: tuple[str]`, `seqlengths: tuple[int\|None]`, `is_circular: tuple[bool\|None]`, `genome: tuple[str\|None]` | Frozen dataclass. `merge(a, b)`: union of seqlevels in R order; conflicting length/circularity/genome → error. |
| `GRanges` | `_seqnames: int32[n]` (codes into `seqinfo.seqnames`), `_ranges: IRanges` (no names/mcols), `_strand: int8[n]` (0 `+`, 1 `-`, 2 `*`, R factor order), `_names`, `_mcols`, `_seqinfo` | Immutable. |
| `Hits` | `query_hits`, `subject_hits: int64` (0-based), `query_length`, `subject_length` | Sorted by query then subject, as R. |
| `Rle` | `values: ndarray`, `lengths: int64` | Internal-grade; only for coverage. |

Python conventions (ADR-0001, ADR-0002):
- Coordinates 1-based like R. **Positional results are 0-based** (Hits, `order`,
  `revmap`, `nearest`/`precede`/`follow`, `select="first"` etc.).
- NA-able integer results: `numpy.ma.MaskedArray[int64]`, masked = NA.
  Never `-1` sentinels (silent last-element indexing).
- API = methods named after R generics in snake_case; kwargs snake_case
  (`min.gapwidth` → `min_gapwidth`, `ignore.strand` → `ignore_strand`).
  Binary ops: `x.intersect(y)`, `x.find_overlaps(y)`. Concat: `IRanges.concat([...])`.
- Setters return copies: `set_names`, `set_mcols`, `set_seqinfo`.
- Subsetting `x[i]`: int, slice, int array, bool mask, list of names.
  Python semantics (0-based, negative from end).
- Warnings R emits (e.g. out-of-bound GRanges) → `warnings.warn(UserWarning)`.

## 3. Public API (method surface)

Shared by IRanges and GRanges unless noted. GRanges adds `ignore_strand` where R has it.

| Group | Methods |
|---|---|
| Construct | `IRanges(start=None, end=None, width=None, names=None, mcols=None)` (any two of three); `GRanges(seqnames, ranges=None, *, start, end, width, strand=None, names=None, mcols=None, seqinfo=None, seqlengths=None)` |
| Access | `start end width names mcols __len__`; GRanges: `seqnames strand ranges seqinfo seqlevels seqlengths is_circular genome` |
| Subset/combine | `__getitem__`, `concat` (GRanges merges seqinfo; mcols columns must match) |
| Intra | `shift(shift)`, `narrow(start, end, width)`, `resize(width, fix="start")`, `flank(width, start=True, both=False)`, `promoters(upstream=2000, downstream=200)`; GRanges: `trim()` |
| Inter | `range()`, `reduce(min_gapwidth=1, drop_empty_ranges=False, with_revmap=False)`, `gaps(start=None, end=None)`, `disjoin()` |
| Set | `union intersect setdiff`, `punion(fill_gap=False)`, `pintersect(resolve_empty="none")`, `psetdiff` |
| Overlaps | `find_overlaps(subject=None, maxgap=-1, minoverlap=0, type="any", select="all")`, `count_overlaps`, `overlaps_any`, `subset_by_overlaps` |
| Nearest | `nearest(subject=None, select="arbitrary")`, `precede`, `follow`, `distance(y)` |
| Coverage | `coverage(shift=0, width=None, weight=1)` → `Rle` (IRanges) / `dict[str, Rle]` in seqlevel order (GRanges) |
| Order | `order(decreasing=False)`, `sort(decreasing=False)`, `rank(ties_method=...)`, `duplicated()`, `unique()` |

Exact defaults and argument sets: copy from R (`args(IRanges::reduce)`,
`showMethods`, `?findOverlaps`). Where R has more arguments than listed,
implement those in TARGET.md; others only if cheap.

## 4. Core algorithms

| Op | Algorithm | Probe R for |
|---|---|---|
| reduce | Stable sort by (group, start, end); sweep, merge while `next.start - cur.end - 1 < min_gapwidth`; revmap = member indices | zero-width ranges absorbed/kept (probe: `reduce(IRanges(c(5,1,20), width=c(0,10,0)))` → `1-10`, `20-19`) |
| gaps | Complement of reduce within [start, end]; default bounds = min start / max end (IRanges), 1 / seqlength (GRanges) | GRanges with seqlengths emits whole-seq gaps for strands with no ranges (probe showed `+` and `-` get `1-20`) |
| disjoin | Breakpoints {start, end+1}; emit segments covered ≥1 | zero-width input |
| setops | `union = reduce(c)`, `intersect`, `setdiff` via reduce + gaps; p* elementwise | `punion` on disjoint (error unless `fill_gap`), `pintersect` empty handling |
| findOverlaps | Sort subject by start; per query `searchsorted` bound on start, filter by end, then type/maxgap/minoverlap predicate. O(n·k) fine; correctness first | zero-width semantics (probe: `5-4` hits `1-10`, `11-10` does not); `maxgap`/`minoverlap` interplay per `type`; strand compatibility (`*` matches all) |
| select | first/last = min/max subject index; arbitrary: match R's actual pick if deterministic, else ADR + membership compare | |
| nearest etc. | distance = gap length, 0 if overlapping/adjacent; nearest prefers overlaps; precede/follow by side; GRanges strand reversal on `-` | tie rules, `*` strand rules, `distance` adjacency |
| coverage | Diff array over [1, width]: +w at start, −w at end+1; cumsum; run-length encode | default width, shift, negative positions, float weight |
| order | `np.lexsort` stable: IRanges (start, width); GRanges (seqlevel, strand, start, width) | `rank` ties methods, `decreasing` stability |
| duplicated/unique | Key tuple as above; first occurrence kept | |

Metadata propagation (names, mcols) per op: probe R for each op and
encode as cases; do not guess. Rule of thumb from R: intra ops, subset,
sort, unique, subsetByOverlaps keep; reduce/gaps/disjoin/range/set ops drop.

## 5. Differential harness

### Case format (YAML, one file per op group)

```yaml
cases:
  - id: reduce-zero-width-absorbed      # unique across all files
    tags: [reduce, zero-width]
    inputs:
      x: {IRanges: {start: [5, 1, 20], width: [0, 10, 0]}}
    call: {fn: reduce, args: [$x], kwargs: {drop.empty.ranges: false}}
  - id: fo-maxgap-within
    inputs:
      q: {IRanges: {start: [1, 10], end: [5, 12], names: [a, b]}}
      s: {GRanges: {seqnames: [chr1], start: [3], width: [4], strand: ["-"],
                    seqinfo: {seqnames: [chr1], seqlengths: [100]}}}
    call: {fn: findOverlaps, args: [$q, $q], kwargs: {type: within, maxgap: 2}}
    compare: exact            # exact (default) | arbitrary
```

- `fn` and kwargs use **R names** (red-team thinks in R). Adapter maps to Python.
- `"$name"` = input reference. Nested call: `{call: {fn: shift, args: [$x, 3]}}`.
- Literals: `IRanges`, `GRanges`, `Seqinfo` objects; `mcols: {col: [..]}`;
  JSON scalars/arrays; `null` = NA. Ambiguous/empty vectors use typed
  form `{int: []}`, `{num: [..]}`, `{chr: [..]}`, `{lgl: [..]}`.
- Escape hatch: `raw: {r: "<expr>", py: "<expr>"}` with inputs bound as
  variables (for `mcols<-`, `seqlevels<-`, etc.). Use sparingly.
- Accepted divergence: `divergence: {adr: "0003", python: <canonical value>}`.
  Harness then asserts Python equals the ADR-sanctioned value (stays a regression test).

### Canonical value (JSON both sides emit)

| Value | Encoding |
|---|---|
| IRanges | `{class: IRanges, start, width, names\|null, mcols\|null}` |
| GRanges | `{class: GRanges, seqnames, start, width, strand, names, mcols, seqinfo: {seqnames, seqlengths, is_circular, genome}}` |
| Hits | `{class: Hits, from, to, nLnode, nRnode}` |
| Rle / RleList | `{class: Rle, values, lengths}` / `{class: RleList, names, elements}` |
| IntegerList (revmap) | `{class: List, elements: [[..],..]}` |
| atomic vector / factor | `{class: vector, type: int\|num\|lgl\|chr, values, names}` |
| error | `{error: "<message>"}` |
| warnings | top-level `warnings: [..]` |

- **Canonical = R's convention: indices 1-based.** R encoder is generic.
  Python adapter declares per-fn `pre`/`post` hooks: `pre` converts R-style
  index args (e.g. `[` with 1-based or negative ints) to Python; `post`
  adds 1 to index-valued results (Hits, order, select≠all, nearest,
  revmap column). Hooks live only in the adapter table; unit-test them.
- R: `simplifyVector=FALSE` on read; `auto_unbox=FALSE`, `na="null"`,
  `digits=NA` on write. Wrap library load in `suppressWarnings`
  (package-built-under warnings appear on stderr).

### Comparison rules
- Structural equality. `int` vs `num` compare by value; `lgl`, `chr` strict.
  Floats: `abs ≤ 1e-9`.
- Errors: pass iff both sides error. Messages not compared.
- Warnings: ignored unless case sets `check_warnings: true`.
- `compare: arbitrary`: Python pick must be a valid R `select="all"` hit (only
  with an ADR justifying it).
- On failure: print R value, Python value, first differing path.

### R execution and cache
- One `Rscript r/run_cases.R in.json out.json` per case file with stale entries.
- `r_cache/<file>.json`: `{case_id: {key, result}}`; key = sha256(canonical
  case JSON + R/IRanges/GenomicRanges versions). Committed.
- Stale key + R available → rerun R. Stale + no R → skip with reason.
  `pytest --refresh-r` forces rerun.
- Version bump (maintenance replay) invalidates all keys → R reruns → `git diff r_cache/` shows behaviour changes.
- Never hand-edit `r_cache/`.

### Maintenance
- New case: append YAML, run `uv run pytest tests/differential -k <id>`, commit case + cache.
- Divergence found: fix Python (preferred) or ADR + `divergence:` entry.
- `tests/differential/README.md` is the red-team guide: case syntax, adapter
  vocabulary, how to extend the adapter for a new `fn`. ≤1 page.

## 6. Staged work (one commit or more per stage)

Every stage: probe R in a scratch Rscript → write cases → implement → `uv run pytest` green → commit (code, cases, cache).

| # | Stage | Done when |
|---|---|---|
| 1 | ADRs 0001/0002 already written. `IRanges` core: construct, accessors, names, mcols, subset, concat; unit tests | unit tests pass |
| 2 | Harness: `run_cases.R`, `encode.R`, codec, adapter, cache, pytest collection, `show.py`, harness README; seed ~15 construction/accessor/subset/concat cases (incl. errors, empty, zero-width, NA names) | cases pass from cold and warm cache |
| 3 | Intra-range IRanges + cases (negative widths error, `fix=center` rounding, `flank(both)`, `promoters` on zero-width) | green |
| 4 | Inter-range: range/reduce/gaps/disjoin + cases (min_gapwidth 0/2, drop_empty, revmap, gaps bounds) | green |
| 5 | Set ops + p* + cases | green |
| 6 | Hits + overlaps (all types, maxgap, minoverlap, select, count/any/subset) + cases | green |
| 7 | nearest/precede/follow/distance + cases (ties, adjacency, empty subject) | green |
| 8 | Rle + coverage + cases | green |
| 9 | order/sort/rank/duplicated/unique + cases | green — **IRanges complete (target ~2h)** |
| 10 | Seqinfo + GRanges core: construct, accessors, subset, concat w/ seqinfo merge, trim, out-of-bound warning + cases | green |
| 11 | GRanges intra (strand-aware resize/flank/promoters; shift not) + cases | green |
| 12 | GRanges inter + set ops, `ignore_strand`, whole-seq gaps + cases | green |
| 13 | GRanges overlaps/nearest/coverage/ordering + cases | green |
| 14 | Sweep: every TARGET op has ≥1 case; resolve divergences; README (design, artifacts, maintenance) | done criteria met |

Stage 1–2 order matters: the harness must exist before stage 3 so every
later op lands with cases. If time runs short, keep breadth: every op
public with ≥1 case before deepening edge cases.

## 7. Risks

| Risk | Mitigation |
|---|---|
| Zero-width semantics differ per op (overlaps, reduce, disjoin, nearest, coverage) | Probe each; ≥2 zero-width cases per op |
| GRanges `gaps`/`reduce` across strands and seqlengths | Probe; cases with unknown and known seqlengths |
| Strand rules for `*` in nearest/precede/follow, flank/resize/promoters on `-` | Probe matrix of strand pairs |
| Adapter index translation hides or invents bugs | Hooks declared per fn, unit-tested; no implicit conversion |
| JSON typing: empty vectors, length-1 unboxing, NA, names NULL vs all-NA | Typed literals; `simplifyVector=FALSE`; `auto_unbox=FALSE` |
| `select="arbitrary"` picks implementation-defined | Try to match; else ADR + `compare: arbitrary` |
| Seqinfo merge order and conflicts on `concat` | Cases with disjoint, overlapping, conflicting seqinfo |
| names/mcols propagation guessed wrong | Canonical value always includes names + mcols |
| Naive O(n·k) overlaps slow | Not a target; note in README; optimize only if a case times out |
| int32 limits, overflow in `shift`/`flank` | Validate like R (error); cases near 2^31 |
| 4h ceiling | Breadth first; IRanges by ~2h checkpoint |

## 8. Artifacts

| Artifact | Owner / when updated |
|---|---|
| `cases/*.yaml` | Implementer seeds; red-team extends. Append-only in spirit; edit only to fix a mistranslated case (note in commit) |
| `r_cache/*.json` | Generated by harness only |
| `adrs/` | Only for accepted divergence or a real design choice. Format: Context, Options, Decision. Linked from cases via `divergence.adr` |
| `NOTES.md` | Append-only, timestamped judgment calls |
| `README.md` | Design summary + how to maintain artifacts (done criterion 3) |
