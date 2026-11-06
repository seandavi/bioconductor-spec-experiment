# Adjudication of non-pass outcomes

Scope: every `mismatch` and `port_error_unexpected` outcome in `results/scored/<port>.jsonl` after the
adapter/harness fixes and rerun (`7aca2ac`): 9,650 cases × 7 ports, **9,378 outcomes in 665 clusters** (34–199
per port). `unsupported` outcomes are not adjudicated. Ports were only read; R 4.6.0 / IRanges 2.47.2 /
GenomicRanges 1.65.1 was re-run for the probes. Independent spot-check: `ledger/reviews/eval-adjudication-spotcheck.md`
(121 clusters, 91.7 % blind agreement on the first adjudication); its split requests are applied below.

## How

| Step | File |
|---|---|
| Join cases, R values, port records, scores | `load.py` |
| Signature per outcome: `port \| op family (* = cross-op argument validation) \| outcome class \| decisive difference \| decisive input feature` | `cluster.py` → `signatures.jsonl` |
| Decisive difference: the first representation normalization under which R and port agree (`nan_as_na`, `num_type`, `names`, `seqlevel_order`, `seqinfo_fields`, `mcols`, `mcol_order`, `row_order`, ...), else `diff` + input feature | `diagnose.py` |
| R self-consistency probes (re-run R on transformed cases with `oracle/oracle.R`) | `probe_coverage_na.py`, `probe_circular.py`, `probe_translate.py`, `probe_seqlevels.py` |
| Is the port also wrong where R is? (`+port_wrong` features) | `cluster.py` `_port_local_na`, `_port_holds_shifted` |
| One rule per signature (regex → class, evidence, rationale) | `verdicts.py` |
| Apply rules, examined cases, per-port totals | `build.py` → `clusters.jsonl`, `summary.json` |
| Q4 tables | `questions.py` → `q4.json` |
| Work lists used while adjudicating | `show.py PORT SUBSTR`, `review.py PORT` |

Reproduce (from `eval/`): `uv run --project . python adjudication/probe_*.py`;
`uv run --project . python adjudication/cluster.py --write`; `uv run --project . python adjudication/build.py --check`
(0 unmatched); `uv run --project . python adjudication/questions.py`. `build/` is gitignored; regenerate with
`just oracle` and `just run PORT`.

`clusters.jsonl`: one line per cluster: `signature, port, n, ops, sources, class, rule, evidence, rationale,
examined (3 ids spread over the cluster, all when n ≤ 3), case_ids`.

### Class rules applied

| Class | Required |
|---|---|
| declared_divergence | a port README / NOTES / ADR / spec states the behaviour differs from R or the feature is not implemented. A spec clause that *claims* R behaviour and is wrong is a port_bug; NOTES marked `UNVERIFIED` are not declarations. |
| r_accident | R evidence that its value is not a contract (equivalent R calls disagree, the value moves under a symmetry, or R returns an object its own `validObject` rejects) **and** the port's value is the principled answer (or an error where the principled answer is an error). |
| adapter_fault | the adapter forwarded a form the port's public API does not take, mistranslated a value, or projected a correct port value into the wrong shape. |
| harness_fault | a defect in the eval harness itself (runner base helpers, projection/comparator contract). |
| port_bug | everything else, including **R accident + the port's own defect on the same case** (declared_divergence if that port behaviour is declared). |

Principled answer, per R-accident kind:

| R accident | Principled answer the port is held to | Check |
|---|---|---|
| NA-weight coverage not local | NA exactly where an NA-weighted range covers, else the sum; R's type, names, lengths | `probe_coverage_na.local_na` on the port value |
| GRanges nearest-family sentinel (coordinates ≤ 0) | R's answer with every range shifted +1000 (R3: its `select="all"` universe) | `probe_translate.jsonl` `r_shifted`, `r_shifted_universe` |
| empty GRanges precede/follow/nearest result class | integer NA vector of query length (IRanges / non-empty answer) | all 64 members checked |
| `resize(ignore.strand=TRUE)` ignores fix vector | R with `ignore.strand=FALSE` (all `+`): starts 33 18 | all 6 members |
| int32 overflow (gaps/setdiff lose the hole, disjoin errors, maxgap non-monotonic, NA coordinates) | int64 arithmetic answer, or an error | each member read |
| `rank()` default errors; `resolve.empty` / `fix` validated lazily | base `rank` average; an error | each member read |

### R probes (evidence for r_accident)

| Probe | Result |
|---|---|
| `probe_coverage_na.py`: coverage with NA weights vs "sum, NA local to the NA-weighted range"; R re-run with ranges reversed | 228 of 372 NA-weight cases put NA where no NA range covers (`coverage(IRanges(c(7,7,6),width=c(1,1,5)), weight=c(2,NA,3))`: positions 8–10 NA, covered only by weight 3); swapping two ranges of a partition turns `NA 3 1 NA 4` into all-NA. |
| `probe_translate.py`: every range +1000 (nearest/precede/follow/distance/overlaps, no circular sequence of known length); now also the shifted `select="all"` universe | 30 of 1,805 cases change, all GRanges nearest-family with coordinates ≤ 0: `heldout/nearest/GRanges/76` returns a chrX query's nearest on **chrY**. Other coord ≤ 0 cases are translation-invariant (port differences there are port bugs). |
| `probe_circular.py`: every `isCircular` TRUE→FALSE | 962 of 1,633 circular-feature cases change; a port value equal to the linear answer = linear semantics. |
| ad hoc `Rscript` (quoted in `verdicts.py`) | `rank(<IRanges>)` default errors; `pintersect(..., resolve.empty=NA\|TRUE\|'a'\|1:2\|integer(0))` errors only once an empty must be resolved; GRanges precede/follow return Hits (nearest: logical) only when a side is empty; `resize(..., fix=c('end','center'), ignore.strand=TRUE)` ignores the fix vector; `resize(<empty>, fix='bogus')` passes; `trim(use.names=FALSE)` keeps names; `findOverlaps(maxgap=2147483645)` hits but 2147483646/7 miss; `promoters/flank` overflow return objects failing `validObject`. |

## Per-port classes (all sources)

| port | adjudicated | clusters | port_bug | declared | r_accident | adapter | harness |
|---|---|---|---|---|---|---|---|
| e0 | 1000 | 71 | 444 | 501 | 55 | 0 | 0 |
| e0b | 1577 | 109 | 1066 | 459 | 52 | 0 | 0 |
| e1 | 490 | 34 | 404 | 83 | 3 | 0 | 0 |
| e2 | 1222 | 86 | 789 | 404 | 29 | 0 | 0 |
| e3 | 1110 | 84 | 657 | 397 | 56 | 0 | 0 |
| e4 | 1308 | 82 | 1105 | 147 | 56 | 0 | 0 |
| biocpy | 2671 | 199 | 2639 | 0 | 29 | 3 | 0 |

## Held-out (headline) and adjusted fidelity

Adjusted = pass / (supported − r_accident − declared − adapter_fault − harness_fault), held-out counts.
"first" = first adjudication (pre-fix results, `58424b1`). Last column: held-out excluding r_accident/adapter/harness
but keeping declared divergences in the denominator.

| port | held-out adjudicated | port_bug | declared | r_accident | adapter | harness | pass/scored (raw) | pass/supported | adj. denom | **adjusted held-out** | first | adjusted, all sources (first) | held-out, declared kept |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| e1 | 82 | 38 | 41 | 3 | 0 | 0 | 2984/3510 (85.0%) | 97.3% | 3022 | **98.7%** | 98.7% | 95.5% (97.9%) | 97.4% |
| e2 | 132 | 94 | 32 | 6 | 0 | 0 | 2891/3510 (82.4%) | 95.6% | 2985 | **96.9%** | 96.8% | 90.7% (90.7%) | 95.8% |
| e4 | 151 | 103 | 21 | 27 | 0 | 0 | 2642/3510 (75.3%) | 94.6% | 2745 | **96.2%** | 96.2% | 86.6% (88.9%) | 95.5% |
| e0 | 170 | 116 | 21 | 33 | 0 | 0 | 2773/3510 (79.0%) | 94.2% | 2889 | **96.0%** | 95.9% | 94.4% (93.9%) | 95.3% |
| e3 | 194 | 157 | 6 | 31 | 0 | 0 | 2571/3510 (73.2%) | 93.0% | 2728 | **94.2%** | 94.2% | 91.4% (91.4%) | 94.0% |
| e0b | 250 | 206 | 15 | 29 | 0 | 0 | 2505/3510 (71.4%) | 90.9% | 2711 | **92.4%** | 92.4% | 86.0% (85.9%) | 91.9% |
| biocpy | 522 | 504 | 0 | 18 | 0 | 0 | 2192/3510 (62.5%) | 80.8% | 2696 | **81.3%** | 80.5% | 66.2% (67.3%) | 81.3% |

Ranking unchanged (e1 > e2 > e4 > e0 > e3 > e0b > biocpy). Held-out moves ≤ 0.8 points; the policy change (R
accident no longer excuses a port defect) bites in the E1 NA-coverage grid, which held-out barely samples: e1
all-source 97.9 → 95.5 %, e4 88.9 → 86.6 %.

## Changes since first adjudication

| Change | Effect |
|---|---|
| Adapter/harness fixes + rerun (`828474e`…`7aca2ac`) | Adjudicated outcomes 10,394 → 9,378. adapter_fault 438 → 3, harness_fault 97 → 0. Forms the ports lack are now `unsupported` (e0/e0b NA dbl/int vectors and gaps bounds, column-name coverage args, `seqinfo=<chr>`, atomic `c()` operand, BiocPy `type="equal"`, scalar-only GRanges coverage/promoters vectors, `incomparables=FALSE`); BiocPy `select=` results project as vectors. |
| R accident + own port defect → port_bug (spot-check policy) | r_accident 1,529 → 280. Of the 1,249 that left: 591 port_bug (e1 226 NaN/double NA coverage; e4 215 NA as 0; BiocPy 144 = 129 NA-coverage crashes / NA as 0, 11 nearest-family values or crashes off R's shifted answer, 1 negative-coordinate follow, 3 int32 (wrap, distance 0); e0b 5, e3 1 nearest-family off R's shifted answer); 226 declared (e2, NOTES.md:27 NA weights refused / not pinned); 432 no longer adjudicated (now `unsupported`: e0 226, e0b 108, BiocPy 98 NA-weight forms). |
| S031 BiocPy `select=` | 224 adapter_fault → 73 pass, 148 port_bug (27 `select_not_applied`: several subjects per query for `first`/`last`/`arbitrary`; 121 in value clusters: circular 64, coord_le0_r_invariant 19, zero_width 18, plain 10, strand 5, names 3, R3 1, int32 1), 3 adapter_fault (below). Matches the spot-check's 76 representation-only / 148 residual split. |
| S041 / S051 / S056 / S050 mixed clusters | Split by the rerun: S056 (27 GRanges NA vector weights) → all `unsupported`. S050 (37 `incomparables=FALSE` crashes) → 18 pass, 19 port_bug value mismatches (BiocPy `unique`/`duplicated`, ordering diff clusters). S051 (36) → 34 GRanges vector-shift `unsupported`, 2 IRanges empty-input crashes port_bug (own rule). S041 (66) → 22 GRanges `unsupported`, 44 IRanges NA-weight crashes port_bug (own rule); 9 of the 44 are R-order-sensitive, but R's value at the given order is the local sum, so the port is held to it. |
| S073 | BiocPy findOverlaps cluster re-signed `hits_not_sorted_by_query` (19): same pairs, Hits sorted by subject; R2 requires query order. Was mis-labelled `nan_as_na` (no-op first ladder step) with a countOverlaps rationale. |
| S116 | BiocPy `mcols(IRanges(names=c('a',NA)))`: port_bug with own rule (NA name stored as 'None'; `use.names=FALSE` is no equivalent call). |
| S014 / S029 | Fixed in the adapters (`9074d01`); no longer adjudicated. |
| int32 splits | `maxgap_int32_max` (R non-monotonic in maxgap) is r_accident for e2/BiocPy and stays declared for e0/e0b/e1/e3/e4; BiocPy `heldout/findOverlaps/GRanges/63` (maxgap 0, R = e0–e4) split out as port_bug. |
| Rules | 9 rules for fixed faults removed; `probe_translate.py` keeps the shifted R3 universe. |
| Side effect of the fixes (not adjudicated) | 144 previously **passing** cases are now `unsupported` (e0 60, BiocPy 34, e0b 20, e2 10, e4 9, e3 6, e1 5), e.g. e0 36 IRanges constructors with NA where R and e0 both error (R1 pass) now "NA in int vector argument". Raw/supported pass counts drop accordingly (e0 held-out pass 2781 → 2773). |

## Ten most consequential port_bug root causes (cases; held-out in brackets)

Clusters with the same port and cause are merged.

| # | n [held-out] | port | clusters | what |
|---|---|---|---|---|
| 1 | 805 [164] | e2 192, e0 141, e3/e4/BiocPy 123, e0b 103 | `*\|r_error_port_value\|flag_validation` | NA/string/number flags (`ignore.strand=NA`, `with.revmap=-1`, `decreasing=1.5`) accepted where R errors "must be TRUE or FALSE"; every port but e1 |
| 2 | 359 [3] | e1 | `coverage\|...(nan_as_na\|num_type)...` (6 clusters) | NA coverage returned as NaN in a double Rle (R NA_integer_/NA_real_), 226 of them where R's NA positions are themselves non-local |
| 3 | 332 [4] | e4 | `coverage\|diff\|r_nonlocal_na+port_wrong` 215, `diff\|na_arg` 117 | NA weights (masked arrays) counted as 0 |
| 4 | 251 [0] | e4 | `overlaps\|equals_R_linear\|circular` | circular overlaps/nearest give R's linear answer; e4 wraps coverage only, declares nothing for overlaps |
| 5 | 233 [0] | e2 | `overlaps\|equals_R_linear\|circular` | no circular handling, undeclared |
| 6 | 210 [3] | BiocPy | coverage NA (4 clusters) | IRanges coverage crashes on NA weights (`int(None)`) or counts them as 0 |
| 7 | 192 [19] | e0b | `subset_concat\|mcols_type_conflict_crash` | `c()` of objects with different mcols columns/types refuses or polars SchemaError; R fills NA and coerces |
| 8 | 173 [0] | BiocPy | `overlaps\|equals_R_linear\|circular` | circular treated as linear |
| 9 | 165 [49] | BiocPy | `coverage\|num_type` | coverage is float64; R integer Rle |
| 10 | 142 [26] | e3 69, e0b 69, e2 4 | `nearest\|ValueError:select must be ...` | `follow(select="last")` / partial `select` rejected (follow takes precede's vocabulary) |

Also large: BiocPy countOverlaps float64 (98 [37]); BiocPy no recycling (94 [11]); BiocPy precede/follow IndexError (93
[13]); e0b p* ops raise for non-identical seqnames (86 [12]).

## Adapter faults remaining (fix not applied)

| Port | Cluster [n] | Cause | Fix |
|---|---|---|---|
| biocpy | `overlaps\|mismatch\|select_result_as_hits` [3]: `e1/overlaps/fo-empty-query-select-last`, `e2/iranges-overlaps/ir.ov.empty.query.first`, `e4/iranges-overlaps/ir-fo-select-empty` | `adapters/biocpy.py:312` `return self._selected_hits(res, n_query) or _Hits(...)`: for an empty query the correct selection `RVector("int", [])` is falsy (`RVector.__len__`, `runner/adapter_base.py:63`), so it falls back to Hits (R: `integer(0)`; port projected Hits 0×n) | `sel = self._selected_hits(res, n_query); return sel if sel is not None else _Hits(res, n_query, n_subject)` |

## Harness faults

None remaining. The first adjudication's harness fault (`RVector.np()` mapping dbl NA to NaN, e0 37 / e0b 60 cases)
is fixed by `828474e` (np() raises on dbl NA; e2/BiocPy opt in with `na_as_nan`). No converter, oracle, projection
or comparator defect was found; R echo self-test 9650/9650 (`results/report.md`).

## Q4: seqlevel order without explicit seqinfo; circular semantics

**Seqlevel order** (`probe_seqlevels.py`, `probe_seqlevels.json`, re-run on the fixed adapters: unchanged):

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
- **BiocPy sorts**: `genomicranges/GenomicRanges.py:186` `SeqInfo(seqnames=sorted(list(set(seqnames))))` —
  lexicographic, so chr10 precedes chr2; it fails both scored first-appearance cases. Its coverage also orders RleList
  names with `np.unique` (sorted) and drops seqlevels without ranges (`e4/granges-coverage/gr-cov-seqlevels-order`).

**Circular semantics** (`questions.py` → `q4.json`; cases where R's answer depends on circularity):

| port | overlaps (match R / linear / other) | nearest | coverage | verdict |
|---|---|---|---|---|
| e1 | 309 / 0 / 3 | 44 / 0 / 0 | 29 / 0 / 0 | implements circular overlaps, nearest and coverage (3 misses = declared adrs/0008-circular-maxgap-hits.md) |
| e4 | 31 / 251 / 10 | 14 / 23 / 7 | 29 / 0 / 0 | coverage only (README.md:67 "coverage (circular wrap)", NOTES.md:16); overlaps/nearest linear, undeclared |
| e0 | 27 / 275 / 10 | 13 / 28 / 3 | 0 / 27 / 0 | none (declared: README.md:102-104 "Not implemented: ... circular-sequence wrap-around") |
| e0b | 31 / 245 / 16 | 14 / 23 / 7 | 1 / 23 / 5 | none (declared: README.md:84,102) |
| e2 | 30 / 233 / 10 | 14 / 27 / 3 | 0 / 29 / 0 | none (undeclared) |
| e3 | 30 / 233 / 10 | 14 / 23 / 7 | 1 / 17 / 1 | none (declared: adr/0014-circular-metadata-only.md) |
| BiocPy | 17 / 175 / 64 | 5 / 15 / 24 | 0 / 0 / 17 | none (linear; "other" adds its float/selection bugs; 12 coverage cases now unsupported) |

Non-zero "match R" for linear ports = circular-feature cases where the circular and linear answers coincide on the
queried pairs. Example: `e1/redteam_circular/rt-circ-origin-coverage` (chrM length 10, range 9–12): R and e1/e4 cover
1–2 and 9–10; e0/e0b/e2/e3 cover only 9–10 (the linear answer).

## Caveats

- Signatures are heuristics; where a mix would have changed the class a feature was added and the cluster split.
  The `+port_wrong` checks are mechanical against the principled answer; the remaining small r_accident clusters
  (empty-result class, resize, rank, resolve.empty, fix, int32) were checked member by member.
- `resolve.empty=1:2|integer(0)` makes e0/e0b/e3/BiocPy raise a numpy "truth value is ambiguous" error rather than a
  validation message; it stays r_accident because the principled (eager) R answer is an error and the comparator
  does not compare error text.
- BiocPy has no divergence docs; its smaller clusters share one evidence/rationale (`verdicts.py` last rule).
- The adapter fixes' broader `unsupported` rules (144 former passes) shrink the supported denominator; that is a
  scoring-scope choice in the adapters, not adjudicated here.

## Post-verification fix (orchestrator)

- Sol's verify note 1: `adapters/biocpy.py` used `selected or Hits`; an empty selected vector is falsy, so 3 empty-query `select=` cases misprojected. Now `is not None`. BiocPy re-scored, `cluster.py --write` + `build.py --check` re-run: 0 unmatched signatures, **adapter_fault = 0 for every port**. BiocPy adjusted held-out unchanged (81.3%).
- Sol's verify note 2 stands: 144 former passes became `unsupported` after the vector-form fixes; part of the adjusted-score movement is scope, not fidelity. Read adjusted fidelity together with raw pass rate.
