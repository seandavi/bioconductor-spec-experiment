# PLAN — e0b-no-r

Plain port of IRanges/GRanges from knowledge alone. No R, no BiocPy, no
harness/goldens/spec. Scope: TARGET.md. Implementer follows stages below,
one commit (or more) per stage.

## 1. Principles

- **Deep kernels, thin classes.** Algorithms operate on numpy arrays
  `(start, end, group)`; classes validate, wrap, carry names/mcols.
  IRanges = one group. GRanges = group code per `(seqlevel, strand)`.
  One kernel serves both.
- **Immutable values.** Every op returns a new object. No setters;
  `with_names()`, `with_mcols()` only.
- **Vectorised numpy.** No per-range Python loops in kernels except over
  groups/seqlevels (few).
- **Python conventions at the API edge, R semantics inside** (table §3).

## 2. Data model

| Type | Fields | Notes |
|---|---|---|
| `IRanges` | `_start: int64[n]`, `_width: int64[n]`, `_names: list[str] \| None`, `_mcols: pl.DataFrame` (height n, may have 0 cols) | `end = start + width - 1`. width ≥ 0. |
| `Seqinfo` | `seqnames: tuple[str]`, `seqlengths: tuple[int\|None]`, `is_circular: tuple[bool\|None]`, `genome: tuple[str\|None]` | `None` = R `NA`. Frozen dataclass. |
| `GRanges` | `_seq: int32[n]` (code into `seqinfo.seqnames`), `_strand: int8[n]` (0 `+`, 1 `-`, 2 `*`), `_ranges: IRanges` (holds names), `_mcols`, `_seqinfo` | Codes give R's factor-level ordering for free. |
| `Hits` | `query_hits: int64[k]`, `subject_hits: int64[k]`, `query_length`, `subject_length` | Sorted by (query, subject). 0-based. |
| `Rle` | `values: ndarray`, `lengths: int64[r]` | Minimal: `__len__`, `to_numpy()`, `__eq__`, `runvalue`/`runlength` aliases. Coverage output only. |

Construction `IRanges(start=None, end=None, width=None, names=None, mcols=None)`:
any two of three; scalars broadcast; all three must agree; width < 0 → `ValueError`.
One constructor only.

`GRanges(seqnames, ranges, strand="*", mcols=None, seqinfo=None)`:
seqnames as list[str] or scalar; strand chars `+ - *`; seqinfo default
= levels in order of first appearance, lengths None. Out-of-bound ranges on
non-circular seqs with known length → `warnings.warn` (R warns, doesn't error).

## 3. API conventions (decide once; log in NOTES)

| Topic | Choice | R |
|---|---|---|
| Coordinates | 1-based closed, as R | same |
| Positional indices (subsetting, Hits, revmap, order, nearest results) | **0-based** | 1-based |
| Missing index (select=first/…, nearest no hit) | `-1` (`NO_HIT` constant; pandas `get_indexer` precedent) | `NA` |
| Negative subscript | Python (from end) | exclusion |
| Arg names | snake_case of R: `min_gapwidth`, `drop_empty_ranges`, `with_revmap`, `ignore_strand`, `maxgap`, `minoverlap`, `select`, `type` | dotted |
| Ops | methods: `x.reduce()`, `x.find_overlaps(y)`, `x.punion(y)` | generics |
| `len(x)` | `__len__` | `length` |
| Integer type | int64, no int32 overflow check | int32 |
| `revmap` | `with_revmap=True` returns `(IRanges, list[np.ndarray])`, 0-based | `revmap` mcols column |

Subsetting `x[i]`: int, slice, int array, bool mask (length must match),
str/list[str] by names (first match). Concatenation:
`IRanges.concat([...])` / `x.concat(y)`; mcols must have identical column
names (R rbind rule) else `ValueError`; names kept if any input has them
(others → `""`).

Metadata survival: intra ops, subsetting, concat, sort keep names+mcols;
inter-range and set ops drop both; overlaps return `Hits`.

## 4. Semantics and algorithms

Definitions: `s,e,w` start/end/width. Empty range: w=0, e=s-1.

### Intra-range (`intra.py`, vectorised, scalar or per-range args)

| Op | Rule |
|---|---|
| `shift(k)` | s+k |
| `narrow(start,end,width)` | relative SEW; negative start/end count from range end (-1 = last pos); result must lie within range else error |
| `resize(width, fix)` | start: s fixed; end: s'=e-w'+1; center: s' = s + (w - w')//2 (`UNVERIFIED` rounding) |
| `flank(width, start=True, both=False)` | start,!both: [s-w, s-1]; end,!both: [e+1, e+w]; both+start: [s-\|w\|, s+\|w\|-1]; both+end: [e-\|w\|+1, e+\|w\|]; negative w, !both: inside flank ([s, s+\|w\|-1] / [e-\|w\|+1, e]) |
| `promoters(upstream=2000, downstream=200)` | IRanges/`+`/`*`: [s-up, s+down-1]; `-`: [e-down+1, e+up]; up/down ≥ 0 |
| GRanges | `-` strand: resize fix start↔end; flank start↔end; promoters as above. shift/narrow not strand-aware |
| `trim()` (GRanges) | non-circular, known length: clip to [1, L]; fully outside → empty range at boundary ([1,0] or [L+1,L]) (`UNVERIFIED`) |

### Inter-range (`inter.py`, kernel on `(start, end, group)`)

- `range`: per group [min s, max e] (includes empty ranges).
- `reduce(min_gapwidth=1, drop_empty_ranges=False, with_revmap=False)`:
  sort by (group, s, e); sweep with running max end `E`; next range joins run
  iff same group and `s - E - 1 < min_gapwidth`. Vectorise via
  `cummax` within groups + `diff` breaks. Empty ranges obey the same rule
  (adjacent/inside → absorbed); drop leftover empties if flag set.
  revmap = original indices per run.
- `gaps(start=None, end=None)`: complement of non-empty reduced ranges in
  [start, end]; defaults min s / max e of x.
- `disjoin()`: R algorithm —
  `starts = unique(s)`, `ends = unique(e)`;
  `adj_start = sorted(unique(starts ∪ ends+1))[:-1]`,
  `adj_end = sorted(unique(ends ∪ starts-1))[1:]`;
  keep pieces overlapping x with minoverlap=1. Per group. Empty ranges add
  breakpoints but never appear.
- GRanges: group = (seq, strand); `ignore_strand=True` → strand forced `*`.
  `gaps` on GRanges: per seqlevel × each of `+ - *`, bounds [1, seqlength];
  seqlength None → max end on that seqlevel (`UNVERIFIED`). So a seqlevel with
  known length and no `-` ranges yields a full-length `-` gap.
- Output order: sorted by (seqlevel code, strand code, start).

### Set ops (`setops.py`)

| Op | Rule (R definitions) |
|---|---|
| `union(y)` | `reduce(concat(x,y), drop_empty_ranges=True)` |
| `setdiff(y)` | with `lo=min starts`, `hi=max ends` of both: `gaps(union(gaps(x,lo,hi), y), lo, hi)` |
| `intersect(y)` | `setdiff(x, gaps(y, lo, hi))` |
| `punion(y, fill_gap=False)` | [min s, max e]; error if gap > 0 and not fill_gap |
| `pintersect(y, resolve_empty="none")` | [max s, min e]; non-overlapping, non-adjacent → error unless `resolve_empty` in `max.start`/`start.x` |
| `psetdiff(y)` | x minus y; error if y strictly inside x (would split) |

GRanges: per (seq, strand), `ignore_strand`, seqinfo merged. Parallel ops
require equal lengths, same seqname; strand mismatch error unless ignore.

### Overlaps (`overlaps.py`)

Predicate for pair (q, s). `gap = max(qs,ss) - min(qe,se) - 1`;
`ov = min(qe,se) - max(qs,ss) + 1`.

| type | hit iff (plus `ov ≥ minoverlap`) |
|---|---|
| any | `gap ≤ maxgap` (default -1 ⇒ must share a position). Empty range at p hits [a,b] iff `a < p ≤ b` when maxgap = -1 (R: "start strictly inside" ⇒ gap -1) (`UNVERIFIED`) |
| start | `\|qs-ss\| ≤ max(maxgap,0)` |
| end | `\|qe-se\| ≤ max(maxgap,0)` |
| equal | both of the above |
| within | `qs ≥ ss and qe ≤ se`; maxgap ignored (`UNVERIFIED`) |

- `type="any"` with both maxgap≠-1 and minoverlap≠0 → `ValueError` (R rule).
- Kernel: sort subject by start; per query, candidate subject window from two
  `searchsorted` calls using `max_subject_width` and slack = max(maxgap,0)+1;
  expand pairs with `np.repeat`; filter by exact predicate. O(n log n + candidates).
  Fine for normal data; pathological for one huge subject (risk §7).
- `select`: all → `Hits`; first/last → min/max subject index; arbitrary →
  same as first (deterministic). Missing → -1.
- GRanges: loop seqlevels present in both; strand compatible iff equal or
  either `*`; `ignore_strand` drops check. `type start/end` on `-`/`-` pairs
  compare the other end (`UNVERIFIED`). Circular wraparound not handled.
- `count_overlaps`, `overlaps_any`, `subset_by_overlaps(invert=False)` derive
  from Hits.

### Nearest (`nearest.py`)

- `distance(y)`: parallel, `max(max(s) - min(e) - 1, 0)`; adjacent = 0.
- `precede(subject, select="first"|"all")`: subjects with `ss > qe`
  (non-overlapping, adjacent allowed), minimal `ss`; tie → lowest index.
- `follow`: subjects with `se < qs`, maximal `se`; tie → lowest index (`UNVERIFIED`).
- `nearest(subject=None, select="arbitrary"|"all")`: overlapping subjects
  win (distance 0); else min distance over precede/follow candidates; tie
  arbitrary → lowest index. `subject=None` → self, excluding self-hits.
- GRanges: same seqlevel, compatible strand. Direction flips on `-`
  (precede = lower coordinates). `*` query treated as `+` (`UNVERIFIED`).
  `distance` across seqlevels/incompatible strand → `-1` (R `NA`), per §3.

### Coverage (`coverage.py`, `rle.py`)

- IRanges `coverage(shift=0, width=None, weight=1)` → `Rle` over positions
  1..width (default max end after shift, 0 if empty). Kernel: `np.add.at` on
  a diff array, cumsum, run-length encode.
- GRanges → `dict[str, Rle]` over all seqlevels; length = seqlength or max end.

### Ordering (`ordering.py`)

- IRanges key: (start, width). GRanges key: (seq code, strand code, start, width).
- `order(decreasing=False)` stable (`np.lexsort`); `sort`; `rank(ties_method=
  "average"|"first"|"min"|"max")` default "average" (`UNVERIFIED` for S4Vectors);
  `duplicated` (first occurrence False); `unique` keeps first, original order.
- No element-wise `__eq__`; `equals(other) -> bool` (starts, widths, names,
  mcols) for tests.

### Seqinfo (`seqinfo.py`)

- `merge(a, b)`: union of seqnames, a's order then new b's; per field `None`
  yields to value; conflicting non-None → `ValueError`. No common seqnames
  (both non-empty) → `warnings.warn`.
- GRanges accessors: `seqlevels`, `seqlengths` (dict), `is_circular`, `genome`.

## 5. Module layout

```
src/ranges_e0b/
  __init__.py    public exports: IRanges, GRanges, Seqinfo, Hits, Rle, NO_HIT
  _sew.py        SEW solving, broadcasting, validation helpers
  iranges.py     IRanges class; methods delegate to kernels
  granges.py     GRanges class; builds group codes, delegates to kernels
  seqinfo.py
  intra.py       shift/narrow/resize/flank/promoters kernels
  inter.py       range/reduce/gaps/disjoin kernels (grouped)
  setops.py
  overlaps.py    Hits + overlap kernel
  nearest.py
  coverage.py
  rle.py
  ordering.py
tests/           test_<module>.py mirroring src
```

## 6. Artifacts and maintenance

| Artifact | Purpose | Maintenance |
|---|---|---|
| `PLAN.md` | this file | update only on redesign; note why in commit |
| `README.md` | design summary, API conventions table (§3), how to run | final stage; keep in sync with §3 |
| `NOTES.md` | append-only, timestamped judgment calls; `UNVERIFIED` tag on R-behaviour guesses | append on every non-obvious semantic choice |
| `tests/` | pytest: hand-derived examples + hypothesis invariants | one file per module |

Test conventions:
- Each example test has a comment with the equivalent R call, e.g.
  `# R: reduce(IRanges(c(1,5,20), c(10,12,25)))`. Lets later evaluation
  replay them in R. Expected values hand-derived; mark doubtful ones
  `# UNVERIFIED`.
- Hypothesis invariants (cheap, catch kernel bugs):
  - reduce output sorted, disjoint, non-adjacent; same covered positions as input.
  - union/intersect/setdiff equal set ops on position sets (brute force).
  - find_overlaps equals brute-force O(nm) predicate for all types/maxgap/minoverlap.
  - coverage sums to total width; disjoin pieces partition reduce.
  - GRanges ops equal IRanges ops applied per (seq, strand).

## 7. Staged work breakdown (≥1 commit each)

| # | Stage | Done when |
|---|---|---|
| 1 | `_sew.py`, `IRanges` core: ctor, accessors, names, mcols, `__len__`, `__repr__`, subsetting, concat, `equals` | tests pass |
| 2 | Intra-range on IRanges | tests incl. negative flank, narrow negatives |
| 3 | Ordering on IRanges | order/sort/rank/duplicated/unique |
| 4 | Inter-range kernels (grouped) + IRanges methods, revmap | invariants + examples |
| 5 | Set ops + p* | brute-force property tests |
| 6 | Overlaps kernel, `Hits`, count/any/subset | brute-force oracle for all types |
| 7 | Nearest family + distance | examples incl. ties, adjacency, self |
| 8 | `Rle` + coverage | examples + sum invariant |
| 9 | `Seqinfo` + merge | conflict/warn tests |
| 10 | `GRanges` core: ctor, accessors, subsetting, concat (seqinfo merge), out-of-bound warning, `trim` | tests |
| 11 | GRanges intra (strand-aware) | `-` strand cases |
| 12 | GRanges ordering | seqlevel order ≠ alpha order case |
| 13 | GRanges inter + set ops, `ignore_strand`, gaps × 3 strands | per-group invariant |
| 14 | GRanges overlaps | strand compatibility matrix |
| 15 | GRanges nearest/distance | strand direction cases |
| 16 | GRanges coverage | dict over seqlevels |
| 17 | README final, `__all__`, done-criteria check (`uv run pytest` green) | TARGET done criteria 1–4 |

Priority if time-boxed (~4 h ceiling): 1–6, 10, 13, 14 first; nearest and
coverage last.

## 8. Risks

| Risk | Mitigation |
|---|---|
| Memory of R edge cases wrong (zero-width overlaps, resize center rounding, nearest ties, `*` direction, gaps with NA seqlength, trim of out-of-range, rank default) | Tag `UNVERIFIED` in tests + NOTES; evaluation adjudicates against R |
| 0-based indices / `-1` sentinel diverge from R's 1-based / NA | Deliberate; documented in README §3 table; projection is mechanical |
| Overlap kernel worst case (one huge subject range → O(nm) candidates) | Acceptable for scope; note in README; NCList is the upgrade path |
| GRanges overlap/strand rules subtle (type start/end on `-`) | Isolated in one function; unit-tested matrix |
| Circular seqs: no wraparound in overlaps/coverage | Out of practical scope; documented |
| polars mcols subsetting/concat quirks (dtype unification) | Strict: identical column names required; use `pl.concat(how="vertical")` |
| Scope creep (setters, Views, lists) | TARGET out-of-scope list; do not add |
