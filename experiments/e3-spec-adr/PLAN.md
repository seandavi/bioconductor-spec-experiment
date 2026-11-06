# PLAN — e3-spec-adr

Spec and ADRs first; code and tests follow them. Reference: R 4.6.0,
IRanges 2.47.2, GenomicRanges 1.65.1, S4Vectors 0.51.6.

## 1. Artifacts

| Artifact | Path | Role | Maintained how |
|---|---|---|---|
| Spec | `spec/NN-area.md` | Normative behavior, language-neutral, clause ids | Edit clause in place; never renumber; retire, don't delete |
| ADRs | `adr/NNNN-slug.md` | Non-obvious choices, linked to clauses | Append-only; supersede with a new ADR |
| Tests | `tests/test_*.py` | Check code against clauses | Each test cites ids via `@pytest.mark.spec(...)` |
| Traceability test | `tests/test_spec_trace.py` | Every live clause cited ≥1 test; every cited id exists | Fails CI on drift |
| NOTES | `NOTES.md` | Judgment calls, timestamped, append-only | Append only |
| README | `README.md` | Design + how to maintain the above | Final stage |

### Spec format

- One file per area (table in §4). Header: area, R versions verified.
- Clause: `**OVL.4** <one or two sentences of semantics>.` then optional
  `Ex:` lines in neutral notation, optional `ADR:` link.
- Notation: ranges as `[start,end]` closed, 1-based; `w` = width;
  zero-width = `[s, s-1]`; `NA` = missing; vectors as `(1, 2, NA)`;
  positions in outputs are 0-based element indices (ADR-0001).
- Describe semantics, not R or Python calls. Argument names in spec
  use the Python spelling (`min_gapwidth`).
- Out of spec = implementation-defined. Do not spec errors/messages
  except where an operation must refuse (e.g. negative width).
- Retired clause: keep id, strike text, add `Retired: <why>`.
- Verify each clause with `Rscript` while writing. Do not commit R
  scripts as goldens; a scratch script under `spec/scratch/` is
  fine if gitignored. Spec examples are the record.

### ADR format

`# ADR-NNNN Title` · Status · Clauses · Context · Options · Decision ·
Consequences. ≤ 1 page. Write one when: Python diverges from R,
R behavior is surprising, or a representation choice leaks into the API.

### Test citation

```python
@pytest.mark.spec("OVL.4", "OVL.5")
def test_zero_width_query_inside_subject() -> None: ...
```

Marker registered in `pyproject.toml`. `test_spec_trace.py` collects
markers (pytest `pytest_collection_modifyitems` hook in `conftest.py`
or AST scan of `tests/`) and regexes `\*\*([A-Z]+\.\d+)\*\*` from
`spec/*.md`.

## 2. Architecture

Deep module, clear seam: all range math lives in pure-numpy kernels
over `(start, width)` int64 arrays. `IRanges` and `GRanges` are thin
classes that validate, dispatch to kernels, and rebuild objects.

### Module layout (`src/ranges_e3/`)

| Module | Contents |
|---|---|
| `_kernels.py` | Pure functions on numpy arrays: intra ops, reduce (+revmap), gaps, disjoin, set ops, overlaps, nearest, coverage, order. No classes. |
| `rle.py` | `Rle` (lengths, values) — minimal, for `coverage` output only. |
| `hits.py` | `Hits` (query_hits, subject_hits, query_length, subject_length). |
| `iranges.py` | `IRanges`. |
| `seqinfo.py` | `Seqinfo`, merge. |
| `granges.py` | `GRanges`; strand/seqname grouping over kernels. |
| `_mcols.py` | mcols helpers (subset, concat, validate height). |
| `__init__.py` | Public exports. |

GRanges may split kernels further (`_kernels/overlaps.py` ...) if
`_kernels.py` passes ~800 lines. Implementer's call; note it.

### Data model

- **IRanges**: `start: int64[n]`, `width: int64[n]` (end derived),
  `names: list[str] | None` (or object array), `mcols: pl.DataFrame | None`.
  Immutable by convention; ops return new objects. Width ≥ 0 enforced.
- **Seqinfo**: `seqnames: list[str]`, `seqlengths: list[int|None]`,
  `is_circular: list[bool|None]`, `genome: list[str|None]`.
- **GRanges**: `seqnames: int32 codes` into `seqinfo.seqnames`,
  `ranges: IRanges` (bare; GRanges owns names and mcols),
  `seqinfo: Seqinfo`, `strand: int8 codes` with order
  `+`=0, `-`=1, `*`=2 (R factor order; drives sorting).
- **Hits**: two int64 arrays, 0-based, sorted by (query, subject).
- **Rle**: `lengths: int64[k]`, `values: int64[k]`. GRanges coverage →
  `dict[str, Rle]` keyed by seqlevel in seqlevel order.

### API conventions (record in ADRs)

| Topic | Recommended decision | ADR |
|---|---|---|
| Index base | Coordinates 1-based closed (as R). Element indices (subsetting, Hits, `nearest`, `order`, revmap) 0-based. | 0001 |
| Missing values | Integer results that may be NA → `numpy.ma.MaskedArray[int64]`. | 0002 |
| mcols | `polars.DataFrame`; empty frame when no columns (polars can't carry height without columns). List columns for `revmap`. | 0003 |
| Naming | R names → snake_case methods (`find_overlaps`, `min_gapwidth`, `ignore_strand`, `with_revmap`). Two-object ops are methods on the first operand. `concat(*objs)` module function. | 0004 |
| Integer range | int64 internally; no int32 overflow emulation. | 0005 |
| Subsetting | `x[int | slice | int array | bool array | str | str list]`; str via names. | 0001 |
| Strand input | Accept `"+"`, `"-"`, `"*"`; reject others. | — |

### Core algorithms (kernels)

- **Intra** (`shift`, `narrow`, `resize`, `flank`, `promoters`):
  vectorized arithmetic. GRanges flips direction for `-` strand
  (resize fix, flank start, promoters upstream). `*` treated as `+`.
- **reduce**: stable sort by (start, end); sweep; merge when
  `next.start <= cur_end + min_gapwidth`; zero-width handling per spec
  (R keeps empties not absorbed unless `drop_empty_ranges`). revmap =
  original indices per output range in R's order (R gave `4,3`: sorted
  order, not input order — spec it).
- **gaps**: complement of reduce within `[start, end]` bounds.
- **disjoin**: breakpoints from all starts and ends+1; keep pieces
  covered by ≥1 range.
- **Set ops**: `union = reduce(c(x,y))`; `intersect`, `setdiff` via
  coverage-like sweep over reduced inputs. `p*` variants elementwise,
  equal length required.
- **Overlaps**: sort subject by start; per query, `searchsorted` for
  candidates with `start <= q_end + maxgap + 1`, then vector mask on
  end/type/minoverlap. O(n·k). Counting may use two searchsorteds.
  Correctness over speed (perf out of scope).
- **Zero-width overlap**: R now reports zero-width `[5,4]` overlapping
  `[1,10]` (verified). Spec the exact rule incl. `minoverlap`.
- **nearest / precede / follow / distance**: sorted starts/ends +
  `searchsorted`; tie rules from R (`select="all"` vs `"arbitrary"`).
- **coverage**: `diff` events at start and end+1, `cumsum`, compress
  to runs. Length = `width` arg, else max end (IRanges) / seqlength
  (GRanges).
- **order**: `np.lexsort` on (end, start) for IRanges; (end, start,
  strand, seqcode) for GRanges. Stable ties. `rank` ties method per R
  default (spec it). `duplicated` first-occurrence.

### GRanges strategy

One helper: group rows by `(seqcode, strandcode)` (or `seqcode` only
if `ignore_strand`), run kernel per group, reassemble.
- Inter-range output order: seqlevel, then strand `+ - *`, then start.
- Overlap compatibility: same seqname; strands compatible if equal or
  either is `*`. `ignore_strand` treats all as `*`.
- `gaps` on GRanges emits all three strands per seqlevel when
  seqlength known (verified); spec unknown-length behavior.
- `trim`: clip to `[1, seqlength]` for non-circular known lengths.
- `concat` merges Seqinfo: union of levels, `x` order then new `y`;
  conflicting lengths/circularity/genome → error.

## 3. Verified R facts (seed the spec)

- `findOverlaps([5,4], [1,10])` → 1 hit, also with `minoverlap=0`.
- `reduce([1,3],[5,4],[10,13],[10,10])` keeps `[5,4]`; revmap `(1),(2),(4,3)` (1-based in R).
- `coverage` of same → runs `3×1, 6×0, 1×2, 3×1`.
- GRanges sort order: seqlevel order (not alphabetical), then `+ - *`, then start, end.
- `findOverlaps` GRanges: `-` query vs `+` subject: no hit; `*` hits both.
- `distance` across strands `+`/`-` → NA; adjacent `[1,3]`,`[4,6]` → 0.
- `resize([1,5], 2, fix=center)` → `[2,3]`; `resize([1,4],3,center)` → `[1,3]`.
- Seqinfo merge: union, `x` levels first; warns on disjoint levels (Python: no warning required).
- `c(GRanges("chr2:1-2"), GRanges("chr3:1-2"))` → seqlevels `chr2, chr3`. Default seqlevel order when built from strings: verify in S4.

## 4. Staged work breakdown

Commit after each row. Message: what + why; end with the
`Co-Authored-By` trailer. Spec amendments discovered during
implementation: separate commit `spec: amend OVL.4 — <why>` + NOTES entry.

| # | Stage | Output | Spec areas / ids |
|---|---|---|---|
| 0 | Scaffold (done by planner) | pyproject, empty pkg, tests | — |
| S1 | Spec: conventions + core | `spec/00-conventions.md`, `10-iranges-core.md` (construction, accessors, names, subsetting, concat, mcols), ADR-0001..0005 | `CONV`, `CORE` |
| S2 | Spec: range ops | `20-intra.md`, `30-inter.md` (range, reduce, gaps, disjoin), `40-setops.md`, `50-coverage.md` | `INTRA`, `INTER`, `SET`, `COV` |
| S3 | Spec: overlaps, nearest, ordering | `60-overlaps.md` (incl. Hits), `70-nearest.md`, `80-order.md` | `OVL`, `NEAR`, `ORD` |
| S4 | Spec: genomic | `90-seqinfo.md`, `95-granges.md` (strand rules per op, `ignore_strand`, trim, sort across seqlevels) + ADRs | `SEQI`, `GR` |
| I1 | IRanges core + trace test | `iranges.py`, `_mcols.py`, `conftest.py`, `test_spec_trace.py` (allow-list of not-yet-implemented areas, shrinks each stage) | CORE |
| I2 | Intra-range | kernels + methods + tests | INTRA |
| I3 | Inter-range + set ops + coverage | kernels, `rle.py`, tests | INTER, SET, COV |
| I4 | Overlaps + Hits | `hits.py`, kernels, tests | OVL |
| I5 | Nearest/precede/follow/distance | kernels, tests | NEAR |
| I6 | Ordering | order/sort/rank/duplicated/unique | ORD |
| I7 | Seqinfo + GRanges core | `seqinfo.py`, `granges.py` construction, accessors, subset, concat+merge, trim | SEQI, GR core |
| I8 | GRanges intra (strand-aware) | | GR intra |
| I9 | GRanges inter/set/coverage | grouping helper | GR inter |
| I10 | GRanges overlaps/nearest/order | | GR rest |
| F | Close-out | README (design + artifact maintenance), trace allow-list empty, `uv run pytest` green | — |

Timebox: S1–S4 ≤ ~75 min total of the 4 h ceiling. Spec only what
an implementer or user needs to predict output; skip print formats,
coercions, error texts.

## 5. Risks

| Risk | Mitigation |
|---|---|
| Spec bloat eats the time budget | Timebox S-stages; one clause per behavior; examples ≤ 2 per clause. |
| Spec drifts from code | Trace test; amend-spec commits; tests cite ids. |
| Zero-width semantics (overlaps, reduce, gaps, disjoin, setops) | Probe R per op in S2/S3; ADR if rule is non-obvious. |
| Nearest ties, `*` strand in nearest/precede/follow | Probe R; spec tie rule per `select`; ADR. |
| `-` strand flips (flank, resize, promoters, precede/follow) | Dedicated GR clauses; table of strand × op. |
| Seqinfo merge conflicts & level order | SEQI clauses; ADR for warning vs silent. |
| revmap / Hits ordering differs from R | Spec ordering explicitly; tests assert order, not sets. |
| 0-based indices vs R's 1-based confuse spec examples | ADR-0001; spec uses 0-based for element indices everywhere. |
| polars zero-column frames lose height | ADR-0003; `_mcols` helper owns the rule. |
| O(n·k) overlaps slow on large inputs | Accept; note in README; out of target scope. |
