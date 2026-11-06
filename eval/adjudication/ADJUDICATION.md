# Adjudication of non-pass outcomes

Scope: every `mismatch` and `port_error_unexpected` outcome in `results/scored/<port>.jsonl` (9,650 cases ×
7 ports; 10,394 outcomes). `unsupported` outcomes are not adjudicated. Ports were only read; R 4.6.0 /
IRanges 2.47.2 / GenomicRanges 1.65.1 was re-run for the probes below.

## How

| Step | File |
|---|---|
| Join cases, R values, port records, scores | `load.py` |
| Signature per outcome: `port \| op family (* = cross-op argument validation) \| outcome class \| decisive difference \| decisive input feature` | `cluster.py` → `signatures.jsonl` |
| Decisive difference: the first representation normalization under which R and port agree (`nan_as_na`, `num_type`, `names`, `seqlevel_order`, `seqinfo_fields`, `mcols`, `mcol_order`, `row_order`, ...), else `diff` + input feature | `diagnose.py` |
| R self-consistency probes (re-run R on transformed cases with `oracle/oracle.R`) | `probe_coverage_na.py`, `probe_circular.py`, `probe_translate.py`, `probe_seqlevels.py` |
| One rule per signature (regex → class, evidence, rationale) | `verdicts.py` |
| Apply rules, examined cases, per-port totals | `build.py` → `clusters.jsonl`, `summary.json` |
| Q4 tables | `questions.py` → `q4.json` |
| Work lists used while adjudicating | `show.py PORT SUBSTR`, `review.py PORT` |

Reproduce (from `eval/`): `python adjudication/probe_*.py`; `uv run --project . python adjudication/cluster.py --write`;
`uv run --project . python adjudication/build.py --check`; `uv run --project . python adjudication/questions.py`.
`build/` (resolved/expected/raw results) is gitignored; regenerate with `just oracle` and `just run PORT`.

`clusters.jsonl`: one line per cluster: `signature, port, n, ops, sources, class, rule, evidence, rationale,
examined (3 ids spread over the cluster, all when n ≤ 3), case_ids`. 700 clusters (38–200 per port); 0 unmatched.

### Class rules applied

| Class | Required |
|---|---|
| declared_divergence | a port README / NOTES / ADR / spec states the behaviour differs from R, or says the feature is not implemented (e.g. "circular wrap-around not implemented"). A spec clause that *claims* R behaviour and is wrong is a port_bug. NOTES marked `UNVERIFIED` are guesses, not declarations. |
| r_accident | R evidence that its value is not a contract: two equivalent R calls disagree, the value moves under a symmetry (translation, range order), or R returns an object its own `validObject` rejects. |
| adapter_fault | the adapter forwarded a form the port's public API does not take (should be `Unsupported`), mistranslated a value, or projected a correct port value into the wrong shape. |
| harness_fault | a defect in the eval harness itself (runner base helpers, projection/comparator contract). |
| port_bug | everything else. |

When R's value is defective the case is r_accident even if the port also has a separate problem there (the
expected value cannot be held against anyone).

### R probes (evidence for r_accident)

| Probe | Result |
|---|---|
| `probe_coverage_na.py`: coverage with NA weights vs "sum, NA local to the NA-weighted range"; R re-run with ranges reversed | 228 of 372 NA-weight cases put NA where no NA range covers (`coverage(IRanges(c(7,7,6),width=c(1,1,5)), weight=c(2,NA,3))`: positions 8–10 NA, covered only by weight 3); swapping two ranges of a partition turns `NA 3 1 NA 4` into all-NA. |
| `probe_translate.py`: every range +1000 (nearest/precede/follow/distance/overlaps, no circular sequence of known length) | 30 of 1,805 cases change, all GRanges nearest-family with coordinates ≤ 0 (sentinel number line): `heldout/nearest/GRanges/76` returns a chrX query's nearest on **chrY**. The other coord ≤ 0 cases are translation-invariant, so port differences there are port bugs. |
| `probe_circular.py`: every `isCircular` TRUE→FALSE | 962 of 1,633 circular-feature cases change; a port value equal to the linear answer = linear semantics. |
| ad hoc `Rscript` (quoted in `verdicts.py`) | `rank(<IRanges>)` default errors (`'arg' must be of length 1`, S4 `.local` passes the choice vector) while `base::rank` defaults to average; `pintersect(x, x, resolve.empty=NA|TRUE|'a'|1:2|integer(0))` returns but errors once an empty must be resolved; GRanges precede/follow return Hits (nearest: logical) only when a side is empty; `resize(..., fix=c('end','center'), ignore.strand=TRUE)` ignores the fix vector (strand-FALSE and IRanges disagree); `resize(<empty>, fix='bogus')` passes; `trim(use.names=FALSE)` keeps names; `promoters/flank` int32 overflow return objects failing `validObject`. |

## Per-port classes (all sources)

| port | adjudicated | clusters | port_bug | declared | r_accident | adapter | harness |
|---|---|---|---|---|---|---|---|
| e0 | 1386 | 86 | 492 | 502 | 281 | 74 | 37 |
| e0b | 1768 | 116 | 1073 | 459 | 165 | 11 | 60 |
| e1 | 496 | 38 | 181 | 83 | 229 | 3 | 0 |
| e2 | 1228 | 89 | 790 | 178 | 255 | 5 | 0 |
| e3 | 1116 | 86 | 659 | 397 | 57 | 3 | 0 |
| e4 | 1316 | 85 | 893 | 147 | 271 | 5 | 0 |
| biocpy | 3084 | 200 | 2476 | 0 | 271 | 337 | 0 |

## Held-out (headline) and adjusted fidelity

Adjusted = pass / (supported − r_accident − declared − adapter_fault − harness_fault), held-out counts.
Last two columns: same formula over all sources; and held-out without excluding declared divergences.

| port | held-out adjudicated | port_bug | declared | r_accident | adapter | harness | pass/scored (raw) | pass/supported | adj. denom | **adjusted held-out** | adjusted, all sources | held-out, declared kept |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| e1 | 87 | 38 | 41 | 5 | 3 | 0 | 2987/3510 (85.1%) | 97.2% | 3025 | **98.7%** | 97.9% | 97.4% |
| e2 | 136 | 95 | 30 | 8 | 3 | 0 | 2896/3510 (82.5%) | 95.5% | 2991 | **96.8%** | 90.7% | 95.9% |
| e4 | 156 | 103 | 21 | 29 | 3 | 0 | 2645/3510 (75.4%) | 94.4% | 2748 | **96.2%** | 88.9% | 95.5% |
| e0 | 187 | 119 | 22 | 35 | 11 | 0 | 2781/3510 (79.2%) | 93.7% | 2900 | **95.9%** | 93.9% | 95.2% |
| e3 | 198 | 158 | 6 | 32 | 2 | 0 | 2573/3510 (73.3%) | 92.8% | 2731 | **94.2%** | 91.4% | 94.0% |
| e0b | 257 | 206 | 15 | 33 | 3 | 0 | 2510/3510 (71.5%) | 90.7% | 2716 | **92.4%** | 85.9% | 91.9% |
| biocpy | 611 | 525 | 0 | 30 | 56 | 0 | 2160/3510 (61.5%) | 78.0% | 2685 | **80.5%** | 67.3% | 80.4% |

Ranking is unchanged by adjudication (e1 > e2 > e4 > e0 > e3 > e0b > biocpy). Held-out has few NA-weight
coverage and circular cases, so r_accident and declared move held-out scores by ≤ 2.3 points; in the E1 suite
they dominate (1,239 NA-weight coverage cases are r_accident across ports; 336/335/289 circular cases are declared
for e0/e0b/e3).

## Ten most consequential port_bug clusters (by cases; held-out share in brackets)

| # | n [held-out] | cluster | what |
|---|---|---|---|
| 1 | 251 [0] | `e4\|overlaps\|mismatch\|equals_R_linear\|circular` | findOverlaps/countOverlaps/overlapsAny on circular sequences give R's *linear* answer; e4 wraps coverage only (README.md:67) and declares nothing for overlaps |
| 2 | 233 [0] | `e2\|overlaps\|mismatch\|equals_R_linear\|circular` | same, e2 has no circular handling and no declaration |
| 3 | 192 [19] | `e0b\|subset_concat\|port_error:call\|mcols_type_conflict_crash` | `c()` of objects with different mcols columns/types: e0b refuses ("columns differ") or polars raises SchemaError; R fills NA and coerces |
| 4 | 192 [33] | `e2\|*\|r_error_port_value\|flag_validation` | NA/string/number flags (`ignore.strand=NA`, `with.revmap=-1`, `decreasing=1.5`) accepted where R errors "must be TRUE or FALSE" |
| 5 | 167 [51] | `biocpy\|coverage\|mismatch\|num_type` | BiocPy coverage is float64; R integer Rle |
| 6 | 142 [0] | `biocpy\|overlaps\|mismatch\|equals_R_linear\|circular` | circular treated as linear |
| 7 | 141 [34] | `e0\|*\|r_error_port_value\|flag_validation` | as #4 for e0 (also 123 [28] each for e3, e4, biocpy; 103 [25] e0b) |
| 8 | 117 [2] | `e4\|coverage\|mismatch\|diff\|na_arg` | NA weights (passed as masked arrays) counted as 0 where R's local answer is NA |
| 9 | 98 [37] | `biocpy\|overlaps\|mismatch\|num_type` | countOverlaps returns float64 |
| 10 | 94 [11] | `biocpy\|*\|port_error:call\|rejects_recycling` | no R recycling/NA-padding of start/width/names/mcols |

Cross-port recurring bugs (not single clusters): flag validation (all ports but e1); GRanges intersect/setdiff not
clipped to [1, seqlength] (every port, 7–15 cases each); countOverlaps names (e0, e0b, e2, e3; declared in e1/e4);
`mcols(<GRanges>)` rownames and `mcols(<IRanges without columns>)` returning a frame instead of NULL; c() formatting
doubles in character columns (`1234567.123456789` vs R's 15 digits); GRanges constructor dropping `mcols(ranges)`;
`follow(select=)` vocabulary (e0b, e3 use precede's `first`; e3's ADR-0010 assumes R's follow takes `first`).

## Adapter faults (fix not applied)

| Port(s) | Cluster(s) [n] | Fault | Fix |
|---|---|---|---|
| biocpy | `overlaps\|mismatch\|select_result_as_hits` [224] | `adapters/biocpy.py:292-296` `_search` projects every BiocFrame as Hits; for `findOverlaps(select="first"/"last"/"arbitrary")` the port's frame already holds one subject per query, R expects an integer vector | when `select != "all"`, project the frame as a length-`n_query` integer vector: subject+1 for queries present, NA otherwise |
| biocpy | `overlaps\|port_error:call\|type_equal_unavailable` [102] | `type="equal"` passed as `query_type` (`_OVERLAP_ARGS`); BiocPy documents only any/start/end/within | raise `Unsupported("type=equal")` |
| biocpy, e0, e0b, e2, e3, e4 | `coverage\|port_error:call\|mcols_column_name_arg` [8, 8, 8, 2, 1, 2] | `coverage(weight="col")` / `shift="col"` forwarded as a string to coverage methods that take numbers only (e0 README API; e0b iranges.py:374; e2 IRanges weight and shifts; e4 shift; BiocPy) | raise `Unsupported("weight/shift as column name")` for those classes, or look the column up (two public calls) |
| all 7 | `construct\|port_error\|AttributeError:'str'/'RVector' ... (seqnames\|_lengths\|get_seqnames)` [2 each] | `GRanges(seqinfo=<character>)` (valid R) forwarded as a str/RVector; every port takes a Seqinfo object | build the port's Seqinfo from the names first, or raise `Unsupported` |
| e0, e0b, e1, e2, e4, biocpy | `subset_concat\|port_error\|c_atomic_operand` [1 each] | `c(<IRanges>, <IRanges>, 1L)` passes the R integer into the port's concat | raise `Unsupported("c: atomic operand")` |
| e0 | `coverage\|mismatch\|(nan_as_na+)num_type\|na_weight_int` [43 + 10] | `adapters/e0.py:64-67` turns an integer vector holding NA into a float64 NaN array (np() would have raised); the double result is the adapter's type change | raise `Unsupported("NA in integer argument")` (as e0b does) |
| e0 | `*\|port_error\|na_arg_crash` [7]; `inter\|mismatch\|diff\|na_arg` [3] | gaps bounds: NA int → NaN (port crashes), NA lgl → None = the port's "argument absent" default (seqlength), not R's "data extent" | raise `Unsupported("NA bound")` |

## Harness faults

| Port(s) | Cluster [n] | Defect | Fix |
|---|---|---|---|
| e0, e0b | `coverage\|mismatch\|nan_as_na\|na_weight_dbl` [37, 60] | `runner/adapter_base.py:82,88` `RVector.np()` maps a **dbl NA to NaN** silently (it raises only for int/lgl/chr; ADAPTERS.md table "dbl NA and NaN both become nan"), and adapters rely on it (`adapters/e0b.py:75-78` expects np() to raise on NA; `e0.py:64-67` mirrors it). spec/PROJECTION.md and the comparator then score the port's NaN against R's NA. The ports' values equal R's answer for the NaN input they received (`coverage(..., weight=c(NaN,3,1,NaN,4))` = NaN 3 1 NaN 4). | make `np()` raise on dbl NA as well (callers wanting NaN opt in, e.g. `np(na_as_nan=True)`), so ports without an NA representation report `unsupported`; or document that NaN carried in from an NA input is not scoreable |

No converter, oracle, projection or comparator defect was found in the examined cases: the R echo self-test is
9650/9650 (`results/report.md:352`), every R value used as evidence was reproduced, and all r_accident evidence is R behaviour, not harness
behaviour.

## Q4: seqlevel order without explicit seqinfo; circular semantics

**Seqlevel order** (`probe_seqlevels.py`, `probe_seqlevels.json`, plus scored cases):

| Input | R | e0 | e0b | e1 | e2 | e3 | e4 | BiocPy |
|---|---|---|---|---|---|---|---|---|
| `GRanges(c("chr2","chr10","chr1","chr2","chrX","chr1"), ...)` | chr2 chr10 chr1 chrX (first appearance) | same | same | same | same | same | same | **chr1 chr10 chr2 chrX** (sorted) |
| `GRanges(c("b","a","c","a"))` | b a c | same | same | same | same | same | same | **a b c** |
| factor seqnames, levels z,a,q, values a z a | z a q (factor levels, unused kept) | a z | a z | a z | **z a q** | a z | a z | a z |
| `seqlengths=c(b=100,a=200)`, seqnames a b | b a | b a | unsupported | **a b** | b a | unsupported | unsupported | unsupported |
| `c(GRanges("chr2"), GRanges(c("chr1","chr10")))` | chr2 chr1 chr10 | same | same | same | same | same | same | same |

- e0, e0b, e1, e2, e3, e4: first appearance, as R (scored: `e2/granges-core/gr.ctor.levels-first-appearance` b c a
  and `e4/granges-construct/gr-con-seqlevels-appearance` chrB chrA chrC pass for all six). Documented: e0 NOTES.md:77,
  e2 NOTES.md:103 and conformance/REPORT.md:37 (`GR.CTOR`), e4 clauses.md:171 (GR-CON-3) and PLAN.md:53.
- Only e2 also follows factor levels (`e2/granges-core/gr.ctor.factor-seqnames-levels`: R z a; e2 z a; others a z).
  Only e0 and e2 follow `names(seqlengths)` order (`e2/granges-core/gr.ctor.seqlengths`, `heldout/GRanges/GRanges/56`:
  R b a / chr2 chr1; e1 a b / chr1 chr2).
- **BiocPy sorts**: `genomicranges/GenomicRanges.py:186` `SeqInfo(seqnames=sorted(list(set(seqnames))))` — lexicographic,
  so chr10 precedes chr2; it fails both scored first-appearance cases (b c a → a b c; chrB chrA chrC → chrA chrB chrC).
  Its coverage also orders RleList names with `np.unique` (sorted) and drops seqlevels without ranges
  (`e4/granges-coverage/gr-cov-seqlevels-order`: R z a, BiocPy a z).

**Circular semantics** (`questions.py` → `q4.json`; cases where R's answer depends on circularity):

| port | overlaps (match R / linear / other) | nearest | coverage | verdict |
|---|---|---|---|---|
| e1 | 309 / 0 / 3 | 44 / 0 / 0 | 29 / 0 / 0 | implements circular overlaps, nearest and coverage (3 misses = declared adrs/0008-circular-maxgap-hits.md) |
| e4 | 31 / 251 / 10 | 14 / 23 / 7 | 29 / 0 / 0 | coverage only (README.md:67 "coverage (circular wrap)", NOTES.md:16); overlaps/nearest linear, undeclared |
| e0 | 27 / 275 / 10 | 13 / 28 / 3 | 0 / 27 / 0 | none (declared: README.md:102-104 "Not implemented: ... circular-sequence wrap-around") |
| e0b | 31 / 245 / 16 | 14 / 23 / 7 | 1 / 23 / 5 | none (declared: README.md:84,102) |
| e2 | 30 / 233 / 10 | 14 / 27 / 3 | 0 / 29 / 0 | none (undeclared) |
| e3 | 30 / 233 / 10 | 14 / 23 / 7 | 1 / 17 / 1 | none (declared: adr/0014-circular-metadata-only.md) |
| BiocPy | 15 / 143 / 105 | 5 / 15 / 24 | 2 / 0 / 27 | none (linear; "other" adds its float/selection bugs) |

"match R" counts that are non-zero for linear ports are circular-feature cases where the circular and linear
answers coincide on the queried pairs. Example: `e1/redteam_circular/rt-circ-origin-coverage` (chrM length 10,
range 9–12): R and e1/e4 cover 1–2 and 9–10; e0/e0b/e2/e3 cover only 9–10 (the linear answer).

## Caveats

- Signatures are heuristics; a few clusters mix root causes that share a class (e.g. e0 `nearest|diff|empty`
  includes an `empty x on a circular sequence` case — both are the same R return-type accident). Where a mix would
  have changed the class, a feature was added and the cluster split (reduce min.gapwidth/drop.empty, distance
  overflow, resolve.empty values R accepts, character c() operands, column-name coverage arguments).
- BiocPy has no divergence docs; its 152 smaller clusters share one evidence/rationale (`verdicts.py` last rule).
- r_accident for NA-weight coverage relies on the local-NA definition; R's own partition-in-order answers (local NA)
  are treated as the contract, so ports that differ there are port bugs (e1: NaN instead of NA; e4/BiocPy: NA as 0).
- Examined ids in `clusters.jsonl` are deterministic picks over the final clusters; the dumps reviewed while
  adjudicating (`show.py`, `review.py`) used the same picker on the then-current clusters.
