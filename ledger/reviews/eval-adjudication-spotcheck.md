# Independent adjudication spot-check

## Final selection (approved revised scope)

- Supervisor explicitly narrowed the original “ALL biocpy” request: largest biocpy clusters covering ≥90% of adjudicated cases, plus 10 random remaining clusters; singleton clusters require their sole member only.
- **121 clusters:** 30 other-port clusters + 91 biocpy clusters. Biocpy largest 81 cover **2,779/3,084 (90.11%)**; 10 random add 21 cases: **2,800/3,084 (90.79%)**. **109 biocpy clusters were not selected.**
- Ordering: largest by `(-n, signature)`; random sampling population sorted lexicographically by signature, excluding largest. Other ports share one `random.Random(7)` stream, in order e0, e0b, e1, e2, e3, e4. Biocpy uses a fresh `Random(7)` for its 10 random choices.
- `L` = largest; `R` = random. This table states every selected signature; seven selected clusters have n=1.

| ID | Port | Pick | n | Signature |
|---|---|---:|---:|---|
| S001 | e0 | L1 | 271 | `e0\|overlaps\|mismatch\|equals_R_linear\|circular` |
| S002 | e0 | L2 | 148 | `e0\|coverage\|mismatch\|nan_as_na+num_type\|r_nonlocal_na` |
| S003 | e0 | L3 | 141 | `e0\|*\|r_error_port_value\|flag_validation\|-` |
| S004 | e0 | R1 | 1 | `e0\|intra\|mismatch\|diff\|strand` |
| S005 | e0 | R2 | 10 | `e0\|construct\|r_error_port_value\|na_input_rejected\|-` |
| S006 | e0b | L1 | 245 | `e0b\|overlaps\|mismatch\|equals_R_linear\|circular` |
| S007 | e0b | L2 | 192 | `e0b\|subset_concat\|port_error:call\|mcols_type_conflict_crash\|-` |
| S008 | e0b | L3 | 106 | `e0b\|coverage\|mismatch\|nan_as_na\|r_nonlocal_na` |
| S009 | e0b | R1 | 24 | `e0b\|inter\|mismatch\|mcols\|-` |
| S010 | e0b | R2 | 15 | `e0b\|ordering\|port_error:call\|ValueError:unknown ties_method 'last'\|-` |
| S011 | e1 | L1 | 118 | `e1\|coverage\|mismatch\|nan_as_na+num_type\|r_nonlocal_na` |
| S012 | e1 | L2 | 108 | `e1\|coverage\|mismatch\|nan_as_na\|r_nonlocal_na` |
| S013 | e1 | L3 | 66 | `e1\|coverage\|mismatch\|nan_as_na\|na_weight_dbl` |
| S014 | e1 | R1 | 2 | `e1\|accessors\|mismatch\|diff\|zero_width` |
| S015 | e1 | R2 | 9 | `e1\|accessors\|mismatch\|names\|-` |
| S016 | e2 | L1 | 233 | `e2\|overlaps\|mismatch\|equals_R_linear\|circular` |
| S017 | e2 | L2 | 192 | `e2\|*\|r_error_port_value\|flag_validation\|-` |
| S018 | e2 | L3 | 118 | `e2\|coverage\|port_error:call\|TypeError:weight must be numeric\|r_nonlocal_na` |
| S019 | e2 | R1 | 1 | `e2\|overlaps\|mismatch\|diff\|int32_edge` |
| S020 | e2 | R2 | 2 | `e2\|construct\|mismatch\|names\|-` |
| S021 | e3 | L1 | 233 | `e3\|overlaps\|mismatch\|equals_R_linear\|circular` |
| S022 | e3 | L2 | 123 | `e3\|*\|r_error_port_value\|flag_validation\|-` |
| S023 | e3 | L3 | 87 | `e3\|*\|port_error:call\|rejects_recycling\|-` |
| S024 | e3 | R1 | 10 | `e3\|nearest\|mismatch\|diff\|empty` |
| S025 | e3 | R2 | 8 | `e3\|setops\|mismatch\|names\|-` |
| S026 | e4 | L1 | 251 | `e4\|overlaps\|mismatch\|equals_R_linear\|circular` |
| S027 | e4 | L2 | 226 | `e4\|coverage\|mismatch\|diff\|r_nonlocal_na` |
| S028 | e4 | L3 | 123 | `e4\|*\|r_error_port_value\|flag_validation\|-` |
| S029 | e4 | R1 | 1 | `e4\|accessors\|mismatch\|diff\|circular` |
| S030 | e4 | R2 | 1 | `e4\|nearest\|r_error_port_value\|na_in_logic\|coord_le0` |
| S031 | biocpy | L1 | 224 | `biocpy\|overlaps\|mismatch\|select_result_as_hits\|-` |
| S032 | biocpy | L2 | 167 | `biocpy\|coverage\|mismatch\|num_type\|-` |
| S033 | biocpy | L3 | 164 | `biocpy\|coverage\|port_error:call\|TypeError:int() argument must be a string, a bytes-like object or a re\|r_nonlocal_na` |
| S034 | biocpy | L4 | 142 | `biocpy\|overlaps\|mismatch\|equals_R_linear\|circular` |
| S035 | biocpy | L5 | 123 | `biocpy\|*\|r_error_port_value\|flag_validation\|-` |
| S036 | biocpy | L6 | 102 | `biocpy\|overlaps\|port_error:call\|type_equal_unavailable\|-` |
| S037 | biocpy | L7 | 98 | `biocpy\|overlaps\|mismatch\|num_type\|-` |
| S038 | biocpy | L8 | 94 | `biocpy\|*\|port_error:call\|rejects_recycling\|-` |
| S039 | biocpy | L9 | 93 | `biocpy\|nearest\|port_error:call\|IndexError:only integers, slices (`:`), ellipsis (`...`), numpy.newaxis\|-` |
| S040 | biocpy | L10 | 84 | `biocpy\|nearest\|mismatch\|diff\|zero_width` |
| S041 | biocpy | L11 | 66 | `biocpy\|coverage\|port_error:call\|TypeError:int() argument must be a string, a bytes-like object or a re\|na_weight_int` |
| S042 | biocpy | L12 | 64 | `biocpy\|seqinfo\|mismatch\|seqinfo_fields\|-` |
| S043 | biocpy | L13 | 63 | `biocpy\|coverage\|mismatch\|diff\|r_nonlocal_na` |
| S044 | biocpy | L14 | 62 | `biocpy\|nearest\|mismatch\|diff\|circular` |
| S045 | biocpy | L15 | 59 | `biocpy\|setops\|mismatch\|seqinfo_fields\|-` |
| S046 | biocpy | L16 | 52 | `biocpy\|setops\|mismatch\|seqinfo\|-` |
| S047 | biocpy | L17 | 45 | `biocpy\|subset_concat\|mismatch\|seqinfo_fields\|-` |
| S048 | biocpy | L18 | 41 | `biocpy\|seqinfo\|r_error_port_value\|seqinfo_incompatible\|-` |
| S049 | biocpy | L19 | 37 | `biocpy\|coverage\|mismatch\|diff\|na_arg` |
| S050 | biocpy | L20 | 37 | `biocpy\|ordering\|port_error:call\|TypeError:argument of type 'bool' is not a container or iterable\|-` |
| S051 | biocpy | L21 | 36 | `biocpy\|coverage\|port_error:call\|Exception:input length N exceeds expected length N\|-` |
| S052 | biocpy | L22 | 33 | `biocpy\|setops\|r_error_port_value\|seqinfo_incompatible\|-` |
| S053 | biocpy | L23 | 32 | `biocpy\|subset_concat\|mismatch\|mcols\|-` |
| S054 | biocpy | L24 | 31 | `biocpy\|overlaps\|mismatch\|diff\|circular` |
| S055 | biocpy | L25 | 28 | `biocpy\|nearest\|mismatch\|diff\|coord_le0_r_invariant` |
| S056 | biocpy | L26 | 27 | `biocpy\|coverage\|port_error:call\|TypeError:int() argument must be a string, a bytes-like object or a re\|na_weight_dbl` |
| S057 | biocpy | L27 | 27 | `biocpy\|setops\|mismatch\|diff\|zero_width` |
| S058 | biocpy | L28 | 26 | `biocpy\|coverage\|mismatch\|diff\|zero_width` |
| S059 | biocpy | L29 | 26 | `biocpy\|nearest\|mismatch\|diff\|strand` |
| S060 | biocpy | L30 | 26 | `biocpy\|setops\|port_error:call\|ValueError:zero-size array to reduction operation minimum which has no \|-` |
| S061 | biocpy | L31 | 24 | `biocpy\|coverage\|mismatch\|diff\|circular` |
| S062 | biocpy | L32 | 24 | `biocpy\|ordering\|mismatch\|diff\|zero_width` |
| S063 | biocpy | L33 | 23 | `biocpy\|inter\|mismatch\|seqinfo_fields\|-` |
| S064 | biocpy | L34 | 23 | `biocpy\|subset_concat\|port_error:call\|ValueError:All objects to combine must have the same number of columns \|-` |
| S065 | biocpy | L35 | 22 | `biocpy\|*\|port_error:call\|rejects_double_arg\|-` |
| S066 | biocpy | L36 | 22 | `biocpy\|nearest\|port_error:call\|ValueError:'query' does not contain the same number of ranges.\|-` |
| S067 | biocpy | L37 | 21 | `biocpy\|ordering\|mismatch\|row_order\|-` |
| S068 | biocpy | L38 | 21 | `biocpy\|overlaps\|r_error_port_value\|seqinfo_incompatible\|-` |
| S069 | biocpy | L39 | 20 | `biocpy\|*\|r_error_port_value\|type_validation\|-` |
| S070 | biocpy | L40 | 20 | `biocpy\|overlaps\|mismatch\|num_type+names\|-` |
| S071 | biocpy | L41 | 20 | `biocpy\|setops\|mismatch\|diff\|circular` |
| S072 | biocpy | L42 | 19 | `biocpy\|construct\|port_error:call\|ValueError:Number of 'seqnames' and 'ranges' do not match!\|-` |
| S073 | biocpy | L43 | 19 | `biocpy\|overlaps\|mismatch\|nan_as_na\|-` |
| S074 | biocpy | L44 | 18 | `biocpy\|accessors\|mismatch\|names\|-` |
| S075 | biocpy | L45 | 17 | `biocpy\|ordering\|mismatch\|diff\|coord_le0` |
| S076 | biocpy | L46 | 16 | `biocpy\|nearest\|mismatch\|diff\|r_translation_dependent` |
| S077 | biocpy | L47 | 16 | `biocpy\|subset_concat\|r_error_port_value\|seqinfo_incompatible\|-` |
| S078 | biocpy | L48 | 15 | `biocpy\|inter\|mismatch\|seqinfo\|-` |
| S079 | biocpy | L49 | 15 | `biocpy\|setops\|mismatch\|diff\|strand` |
| S080 | biocpy | L50 | 14 | `biocpy\|intra\|r_error_port_value\|other:solving row N: 'allow.nonnarrowing' is FALSE and the supplie\|-` |
| S081 | biocpy | L51 | 14 | `biocpy\|nearest\|mismatch\|equals_R_linear\|circular` |
| S082 | biocpy | L52 | 13 | `biocpy\|construct\|mismatch\|names\|-` |
| S083 | biocpy | L53 | 13 | `biocpy\|ordering\|mismatch\|diff\|strand` |
| S084 | biocpy | L54 | 12 | `biocpy\|construct\|r_error_port_value\|other:'start' and 'width' must be numeric vectors\|-` |
| S085 | biocpy | L55 | 12 | `biocpy\|nearest\|r_error_port_value\|other:'arg' must be NULL or a character vector\|-` |
| S086 | biocpy | L56 | 11 | `biocpy\|coverage\|mismatch\|diff\|strand` |
| S087 | biocpy | L57 | 11 | `biocpy\|metadata\|mismatch\|names\|-` |
| S088 | biocpy | L58 | 11 | `biocpy\|nearest\|mismatch\|diff\|empty` |
| S089 | biocpy | L59 | 10 | `biocpy\|*\|r_error_port_value\|scalar_validation\|-` |
| S090 | biocpy | L60 | 10 | `biocpy\|coverage\|mismatch\|diff\|coord_le0` |
| S091 | biocpy | L61 | 10 | `biocpy\|setops\|mismatch\|diff\|coord_le0` |
| S092 | biocpy | L62 | 9 | `biocpy\|construct\|port_error:call\|ValueError:Length of 'strand' does not match the number of ranges.\|-` |
| S093 | biocpy | L63 | 9 | `biocpy\|intra\|port_error:call\|TypeError:only integer scalar arrays can be converted to a scalar inde\|-` |
| S094 | biocpy | L64 | 9 | `biocpy\|overlaps\|mismatch\|diff\|coord_le0_r_invariant` |
| S095 | biocpy | L65 | 9 | `biocpy\|overlaps\|mismatch\|diff\|zero_width` |
| S096 | biocpy | L66 | 9 | `biocpy\|setops\|port_error:call\|KeyError:'chrN'\|-` |
| S097 | biocpy | L67 | 8 | `biocpy\|coverage\|port_error:call\|mcols_column_name_arg\|-` |
| S098 | biocpy | L68 | 8 | `biocpy\|intra\|r_error_port_value\|other:solving row N: either the supplied start or the supplied end\|-` |
| S099 | biocpy | L69 | 8 | `biocpy\|metadata\|port_error:call\|ValueError:Number of 'seqnames' and 'ranges' do not match!\|-` |
| S100 | biocpy | L70 | 8 | `biocpy\|overlaps\|mismatch\|diff\|plain` |
| S101 | biocpy | L71 | 8 | `biocpy\|setops\|mismatch\|mcols\|-` |
| S102 | biocpy | L72 | 7 | `biocpy\|construct\|mismatch\|mcols\|-` |
| S103 | biocpy | L73 | 7 | `biocpy\|construct\|mismatch\|seqlevel_order\|-` |
| S104 | biocpy | L74 | 7 | `biocpy\|coverage\|port_error:call\|RuntimeError:Unable to cast Python instance of type <class 'list'> to C++\|-` |
| S105 | biocpy | L75 | 7 | `biocpy\|inter\|r_error_port_value\|seqinfo_incompatible\|-` |
| S106 | biocpy | L76 | 7 | `biocpy\|nearest\|port_error:call\|IndexError:only integers, slices (`:`), ellipsis (`...`), numpy.newaxis\|r_translation_dependent` |
| S107 | biocpy | L77 | 7 | `biocpy\|nearest\|r_error_port_value\|seqinfo_incompatible\|-` |
| S108 | biocpy | L78 | 7 | `biocpy\|ordering\|mismatch\|diff\|plain` |
| S109 | biocpy | L79 | 7 | `biocpy\|seqinfo\|mismatch\|diff\|coord_le0` |
| S110 | biocpy | L80 | 7 | `biocpy\|setops\|mismatch\|diff\|out_of_bounds` |
| S111 | biocpy | L81 | 6 | `biocpy\|coverage\|mismatch\|diff\|multi_seqlevel` |
| S112 | biocpy | R1 | 1 | `biocpy\|intra\|mismatch\|diff\|strand` |
| S113 | biocpy | R2 | 2 | `biocpy\|coverage\|port_error:call\|RuntimeError:Unable to cast Python instance of type <class 'float'> to C+\|-` |
| S114 | biocpy | R3 | 2 | `biocpy\|intra\|r_error_port_value\|other:'upstream' must be a vector of integers\|-` |
| S115 | biocpy | R4 | 1 | `biocpy\|overlaps\|r_error_port_value\|other:when 'type' is "any", at least one of 'maxgap' and 'minoverl\|-` |
| S116 | biocpy | R5 | 1 | `biocpy\|accessors\|r_error_port_value\|na_input_rejected\|-` |
| S117 | biocpy | R6 | 4 | `biocpy\|construct\|port_error:call\|TypeError:'strand' must be either a numpy vector, a list of integers o\|-` |
| S118 | biocpy | R7 | 3 | `biocpy\|setops\|port_error:call\|KeyError:'chrY'\|-` |
| S119 | biocpy | R8 | 2 | `biocpy\|ordering\|mismatch\|diff\|circular` |
| S120 | biocpy | R9 | 3 | `biocpy\|construct\|port_error:call\|ValueError:end position should fit in a N-bit signed integer\|-` |
| S121 | biocpy | R10 | 2 | `biocpy\|intra\|r_error_port_value\|other:'downstream' must be a vector of integers\|-` |

## Method / independence

- Read the four requested eval documents, adapter contract, port READMEs/NOTES/public source. Initially extracted only signature, port, n, case_ids from adjudication JSONL; class vocabulary and class-rule definitions were read, but **individual cluster class/rationale/evidence remained hidden until blind labels were saved**. The class-rules document also exposed its general R-probe summary; the symmetry probes below were independently written and executed, not copied from those scripts.
- Per selected cluster: freshly ran lexicographically first and last member (sole member for n=1) through live R and the frozen port adapter. **235 case-port runs; 201 distinct R cases.** All 201 fresh R values and all 235 port outcomes reproduced stored records exactly, excluding timing/warnings.
- R receives original literals from `eval/cases/{e1,e2,e4,heldout}.jsonl`; Python receives `eval/build/resolved.jsonl`. An initial attempt to send resolved projections to R produced `r_input_error`; those records were discarded and replaced by the successful literal-input rerun.
- Blind decisions frozen at 2026-11-06 15:21:07 UTC, `/tmp/spotcheck/own-frozen.json`, SHA256 `ef030faf0fd181262ea6ce9171907e266b6ffee44e425c411085d97b60e962e8`. Only then opened cluster verdicts/rationales.
- Additional verification: coverage auto/sort/hash/naive, range+weight reversal, nearest-family +1000 translation, resize strand-flag equivalence, maxgap monotonicity, empty-result classes, and native API signatures. Expanded biocpy reruns covered all members of S031/S041/S050/S051/S056/S093/S104 (406 cases). Including initial probes: **627 unique case-port pairs / 593 unique R cases**.
- Runs used existing venv interpreters directly, `-B` / `PYTHONDONTWRITEBYTECODE=1`, avoiding uv sync or port changes. No repository files changed, no commits, no staged files; all six frozen worktrees remained clean at their listed SHAs.

## Blind agreement

Agreement means every sampled member’s independently recorded class equals the cluster verdict. A mixed sample is a cluster disagreement. Class denominator is the **adjudicator’s** class; unweighted clusters/cases, not projected population accuracy.

| Class | Clusters agreeing / sampled | Rate | Case-port labels agreeing / sampled | Rate |
|---|---:|---:|---:|---:|
| port_bug | 88/97 | 90.72% | 176/191 | 92.15% |
| declared_divergence | 5/5 | 100% | 10/10 | 100% |
| r_accident | 16/16 | 100% | 28/28 | 100% |
| adapter_fault | 2/3 | 66.67% | 4/6 | 66.67% |
| harness_fault | 0 sampled | n/a | 0 sampled | n/a |
| **Overall** | **111/121** | **91.74%** | **218/235** | **92.77%** |

| Port | Cluster agreement |
|---|---:|
| e0 | 5/5 |
| e0b | 5/5 |
| e1 | 4/5 |
| e2 | 5/5 |
| e3 | 5/5 |
| e4 | 4/5 |
| biocpy | 83/91 |

- These are **before-rationale** rates. Two initial positions were subsequently revised below; the blind rate is deliberately not retroactively improved.

## Every blind disagreement

PB = port_bug; DD = declared_divergence; RA = r_accident; AF = adapter_fault. Case IDs and original per-member labels are in the audit register below.

| Cluster | Blind → adjudicator | Fresh evidence / comparison / disposition |
|---|---|---|
| S014, e1 (n=2) | AF → PB | Both `heldout/mcols/GRanges/{46,53}`: R returns 3×0 / 5×0 DataFrames; adapter returns NULL. `ranges_e1/granges.py:206–207` uses None for absent columns; `adapters/e1.py:345–346` forwards it without owner length/names. **Retain representation-boundary objection:** documented storage sentinel can be projected with owner context; adjudicator instead treats missing GRanges-specific accessor shape as semantic. This is a contested boundary, not proof that native accessor fidelity is correct. |
| S029, e4 (n=1) | AF → PB | `e2/granges-core/gr.mcols.none`: R 6×0 DataFrame; port projection 0×0. `ranges_e4/granges.py:121–122` returns the stored zero-column polars frame; `adapters/e4.py:272–273` loses owner length. `adapters/e4.md` explicitly documents the representation limit. **Same boundary objection as S014**; zero-column storage cannot itself preserve height. |
| S031, biocpy (n=224) | PB → AF | Initial two cases only exposed Hits-vs-vector shape. `_search` (`biocpy.py:292–296`) does misproject selected hits, and the correct selected frame can be padded/projected without rerunning search. **Accept AF for representation-only members. However, full-cluster rerun shows 148/224 remain wrong after that representation repair** (details below). Do not excuse all 224 as adapter-only. |
| S041, biocpy (n=66) | RA / AF → PB | `covna0-partition-int`: R runs `NA,3,1,NA,4`, lengths `1,3,3,3,2`; reverse ranges **and weights** → all-NA length12. Port raises int(None) TypeError. This satisfies the adjudicator’s own range-order RA test, although the original expected value itself is sensible. Other sampled member `heldout/coverage/GRanges/73` forwards a two-element weight list to scalar-only GRanges.coverage. **Mixed**: expanded rerun finds 14 order-sensitive members; among 52 invariant members, 35 IRanges / 17 GRanges. Under stated priority: RA / PB / AF respectively. Genuine NA handling deficiency remains; see caution below. |
| S050, biocpy (n=37) | AF → PB | Both cases explicitly pass R `incomparables=FALSE`; adapter `_val` sends Python bool to `biocutils.duplicated/unique`, whose public API requires a membership-capable container (`biocutils/duplicated.py:12–24`). That explains the exact bool-not-iterable error. Calling the native generic with `incomparables=[]` removes this crash. **Retain AF for forwarding**, but this is not exoneration: `unique` still keeps a repeated `(start=0,width=6)` in the second case after that repair. Separate port defect exists. |
| S051, biocpy (n=36) | AF / PB → PB | `gr-coverage-shift-per-range`: R accepts four per-range shifts; native `GenomicRanges.coverage(shift:int, width:int\|None, weight:int)` has no vector-shift API. `heldout/coverage/IRanges/41`: empty IRanges with scalar shift4 is valid for documented array-valued IRanges API after thin conversion, yet native normalization raises “input length 1 exceeds expected length 0”. **Retain mixed classification**: all-member rerun gives 34 GRanges form-limit cases / 2 genuine IRanges empty-input failures. |
| S056, biocpy (n=27) | AF → PB | Both sampled cases are GRanges with multi-element double weight lists containing NA; both crash int(None). Scalar-only `GenomicRanges.coverage` was given an unsupported vector form; adapter already notes this scalar/array distinction (`biocpy.py:569–573`). **Retain AF for these members**, not a recovered pass. Expanded reversal additionally finds 2 RA members; remaining25 are GRanges API-form cases under the written AF rule. |
| S093, biocpy (n=9) | AF → PB | `gr-promoters-vector` and `heldout/promoters/GRanges/50`: R accepts per-range upstream vectors; native `GenomicRanges.promoters(upstream:int, downstream:int)` documents scalars (`GenomicRanges.py:1326–1338`), adapter forwards lists, native raises integer-scalar-index TypeError. **Retain AF under the explicit “form the public API does not take → Unsupported” rule**; all9 members use GRanges. If unsupported vector features are instead intended PB, amend that rule consistently. |
| S104, biocpy (n=7) | AF → PB | Both width-vector cases pass `[5,6,7]` to scalar `GenomicRanges.coverage(width:int\|None)`; native raises list-to-C++ RuntimeError. **Retain AF under the same rule**; all7 members are GRanges. R’s vector-width feature is absent, not merely a broken accepted native scalar call. |
| S116, biocpy (n=1) | RA → PB | R accepts IRanges names `a,NA`, but default `mcols` errors because DataFrame rownames cannot contain NA. `mcols(x,use.names=FALSE)` succeeds; that flag changes requested semantics and is **not an equivalent-call proof**. Moreover BiocPy already converted NA name to literal `"None"`. **Withdraw my initial RA label; accept PB.** A surprising R error alone does not satisfy RA. |

## Adapter over-excusal: S031 must not be treated as wholly repaired

- Full cluster: reran all224 in R and BiocPy. In /tmp only, converted selected native hit pairs to a length-nquery 1-based vector, padding absent queries with NA. No recomputation, sorting, or altered port method.
- First/last: **36/148** now exact; **112/148** still differ. Arbitrary: fresh `select="all"` universes and R3’s union with R’s own pick give **40/76** admissible; **36/76** still fail. Total: **76 representation-only / 148 residual semantic mismatches**. The counts describe observed failures after normalization, not 148 newly proven homogeneous PB verdicts.
- Non-circular control `e1/granges_ops/gr-fo-select-first`: R `[1,2,NA,5]`; normalized native `[1,NA,NA,5]`. Query2 is chr1:10–17 on `-`; subject2 is chr1:12–14 on `*`. That overlap exists; padding cannot recover a hit the port did not emit.
- Representation-only control `gr-fo-empty-subject-select-first`: native empty hit frame with nquery1 converts exactly to R `[NA]`.
- **Split representation-only from co-fault cases**, or retain separate adapter/port fault flags. A single AF bucket excluded from adjusted fidelity conceals a genuine native failure demonstrated above.

## Other mixed or misleading signatures

| Signature / cluster | Finding |
|---|---|
| S041 `coverage / int(None) / na_weight_int` | One exception signature merges IRanges vs scalar-only GRanges APIs and 14 R-order-sensitive cases. These need distinct primary classes under the written rules; the expected local value is not itself wrong in every RA-labelled symmetry test. |
| S051 `coverage / input length N exceeds expected length N` | Different cause: vector-to-scalar API forwarding vs empty-array native bug (34/2). |
| S056 `coverage / int(None) / na_weight_dbl` | All27 GRanges API-form cases; 2 additionally fail the R range-order symmetry, giving RA precedence under current rules. |
| S050 `ordering / bool not iterable` | AF explains the actual scored crash, but fixing it uncovers independent native duplicate-detection failure. Do not convert AF into a claim of correct port semantics. |
| S073 `overlaps / mismatch / nan_as_na` (n=19) | Class PB agrees, but rationale/evidence says “countOverlaps type / float64”. Cluster ops are **19 findOverlaps**, and both fresh samples have **unsorted query hit indices**, not a count type problem. E.g. heldout findOverlaps/IRanges/58: R from `[2,3,4,4,5,5]`; port `[2,3,4,5,4,5]`. R2 requires nondecreasing query order. Rewrite explanation/signature; same class, different cause. |

## r_accident checks / policy cautions

| Sampled RA mechanism | Independent check | Assessment |
|---|---|---|
| Nonlocal NA coverage: S002/S008/S011/S012/S018/S027/S033/S043 | On each sampled endpoint, fresh R `method="naive"` differs from auto; numeric/NA weights are unchanged. In `heldout/coverage/IRanges/14`, shifted NA range is positions11–13, but auto stays NA through28; naive is local. | Evidence meets equivalent-algorithm criterion. Not every range reversal changes this output; no blanket “all NA coverage is RA” claim. |
| Resize: S004/S112 | Same plus-strand GRanges, width3, fix `[end,center]`: `ignore.strand=TRUE` yields starts `[33,23]`; FALSE yields `[33,18]`, the port answer. Strands all `+`, so flag should not change anchoring. | RA supported. |
| Overflow: S019 | R query1:1 vs subject2147483642:2147483642: maxgap2147483645 hits; increasing gap to2147483646 or2147483647 loses hit. E2 finds hit. | Non-monotonic threshold demonstrates arithmetic defect, not merely “large inputs look odd”. |
| Empty shape: S024/S088 | Fresh empty-subject nearest gives IRanges integer vs GRanges logical; precede/follow give IRanges integer vs GRanges Hits. | Consistency criterion supports RA; distinguish this from a substantive missing nearest hit. |
| Sentinel: S030/S076/S106 | +1000 to both operands changes fresh R follow/nearest/precede results or removes logical-NA error. S055/S094 controls stay unchanged under same translation. | RA supported only for demonstrated affected cases; coordinates≤0 alone are insufficient. |

- **No sampled adjudicator-RA cluster was rejected**, but this does not validate every member of those16 clusters. Conversely my initial S116 RA was overreach and is withdrawn.
- Current RA precedence explicitly ignores separate port defects whenever R fails a consistency test. S033, for example, still crashes on NA; S002/S011 reproduce bad NA positions but also use the wrong NA/type encoding. Excluding the whole case flatters a port’s adjusted score even when the R-defective leaf was not the source of its scored failure.
- S041 is the reverse caution: a correct local expected value becomes order-sensitive only through an equivalent permutation that exposes R’s bug. The written rule supports RA, but blindly excusing a native NA crash overstates what that evidence establishes.
- AF should not become a pass. For documented missing forms, Unsupported still belongs in the raw denominator. For S031’s co-faults, do not exclude all native errors simply because a projection problem also exists.

## Audit register: blind labels and fresh cases

- A/B are the first/last cluster member; n=1 has A only. Entries below were run before individual verdict rationales. Slash-separated blind classes correspond to A/B; the adjudicator applies its one class to every member.

| Cluster | A | B | Blind | Adjudicator |
|---|---|---|---|---|
| S001 | `e1/divergences/div-circ-end-maxgap-context` | `e1/redteam_circular/rt-circ-zero-interior` | DD | DD |
| S002 | `e1/grid_coverage_na/covna0-partition-int-gr2` | `heldout/coverage/IRanges/14` | RA | RA |
| S003 | `e1/grid_bad_args/bad-gr-disjoin-revmap-NA` | `heldout/union/GRanges/30` | PB | PB |
| S004 | `heldout/resize/GRanges/35` | — n=1 | RA | RA |
| S005 | `e1/granges_core/gr-seqnames-na-error` | `e4/granges-edge/gr-edge-seqname-empty-string` | PB | PB |
| S006 | `e1/divergences/div-circ-end-maxgap-context` | `e1/redteam_circular/rt-circ-zero-interior` | DD | DD |
| S007 | `e1/construct/ir-concat-mcols-different-columns` | `heldout/c/IRanges/8` | PB | PB |
| S008 | `e1/grid_coverage_na/covna1-partition-num-gr2` | `heldout/coverage/IRanges/44` | RA | RA |
| S009 | `e1/grid_granges/grid-gr-r1-reduce-g1-drop-ign` | `e4/iranges-inter/ir-reduce-revmap-drop-empty` | PB | PB |
| S010 | `e1/grid_granges/grid-gr-rank-last-a` | `heldout/rank/IRanges/6` | PB | PB |
| S011 | `e1/grid_coverage_na/covna0-partition-int-gr2` | `heldout/coverage/IRanges/14` | RA | RA |
| S012 | `e1/grid_coverage_na/covna1-partition-num-gr2` | `heldout/coverage/IRanges/44` | RA | RA |
| S013 | `e1/grid_coverage_na/covna1-partition-num` | `e1/grid_coverage_na/covna9-partition-num-shift` | PB | PB |
| S014 | `heldout/mcols/GRanges/46` | `heldout/mcols/GRanges/53` | AF | PB |
| S015 | `e1/granges_core/gr-mcols-accessor` | `heldout/mcols/IRanges/21` | DD | DD |
| S016 | `e1/divergences/div-circ-end-maxgap-context` | `e1/redteam_circular/rt-circ-zero-interior` | PB | PB |
| S017 | `e1/grid_bad_args/bad-gr-disjoin-revmap-NA` | `heldout/union/GRanges/30` | PB | PB |
| S018 | `e1/grid_coverage_na/covna0-partition-int-gr2` | `heldout/coverage/IRanges/14` | RA | RA |
| S019 | `e4/iranges-edge/ir-edge-fo-maxgap-int-max` | — n=1 | RA | RA |
| S020 | `e1/grid_construct/ctor-ir-names-dblwhole` | `e1/grid_construct/ctor-ir-names-lgl` | PB | PB |
| S021 | `e1/divergences/div-circ-end-maxgap-context` | `e1/redteam_circular/rt-circ-zero-interior` | DD | DD |
| S022 | `e1/grid_bad_args/bad-gr-disjoin-revmap-NA` | `heldout/union/GRanges/30` | PB | PB |
| S023 | `e1/coverage/cov-shift-recycle-warns` | `heldout/resize/GRanges/43` | DD | DD |
| S024 | `e1/granges_ops/gr-nearest-empty-subject` | `heldout/precede/GRanges/31` | RA | RA |
| S025 | `heldout/pintersect/IRanges/18` | `heldout/punion/IRanges/3` | PB | PB |
| S026 | `e1/divergences/div-circ-end-maxgap-context` | `e1/redteam_circular/rt-circ-zero-interior` | PB | PB |
| S027 | `e1/grid_coverage_na/covna0-partition-int-gr2` | `heldout/coverage/IRanges/44` | RA | RA |
| S028 | `e1/grid_bad_args/bad-gr-disjoin-revmap-NA` | `heldout/union/GRanges/30` | PB | PB |
| S029 | `e2/granges-core/gr.mcols.none` | — n=1 | AF | PB |
| S030 | `e1/redteam_edges/rt-edge-gr-negative-follow` | — n=1 | RA | RA |
| S031 | `e1/granges_ops/gr-fo-empty-subject-select-first` | `heldout/findOverlaps/IRanges/9` | PB | AF |
| S032 | `e1/coverage/cov-adjacent-ranges-merge-runs` | `heldout/coverage/IRanges/8` | PB | PB |
| S033 | `e1/grid_coverage_na/covna0-partition-int-gr2` | `heldout/coverage/IRanges/14` | RA | RA |
| S034 | `e1/divergences/div-circ-end-maxgap-context` | `e1/redteam_circular/rt-circ-zero-interior` | PB | PB |
| S035 | `e1/grid_bad_args/bad-gr-disjoin-revmap-NA` | `heldout/union/GRanges/30` | PB | PB |
| S036 | `e1/grid_circular/circ-equal0-fo-equal` | `heldout/subsetByOverlaps/IRanges/27` | AF | AF |
| S037 | `e1/grid_bad_args/bad-gr-count-maxgap-neg` | `heldout/countOverlaps/IRanges/9` | PB | PB |
| S038 | `e1/construct/ir-names-too-short-padded` | `heldout/set_names/IRanges/29` | PB | PB |
| S039 | `e1/granges_ops/gr-follow` | `heldout/precede/GRanges/59` | PB | PB |
| S040 | `e1/grid_granges/grid-gr-distance-equal` | `heldout/nearest/GRanges/69` | PB | PB |
| S041 | `e1/grid_coverage_na/covna0-partition-int` | `heldout/coverage/GRanges/73` | RA/AF | PB |
| S042 | `e1/granges_core/seqinfo-merge-basic` | `heldout/merge/Seqinfo/9` | PB | PB |
| S043 | `e1/grid_coverage_na/covna11-partition-num-shift` | `heldout/coverage/IRanges/6` | RA | RA |
| S044 | `e1/grid_circular/circ-equal0-nearest-all` | `e1/grid_seqinfo/si6-distance` | PB | PB |
| S045 | `e1/grid_seqinfo/si0-setdiff` | `heldout/union/GRanges/48` | PB | PB |
| S046 | `e1/granges_inter/gr-intersect-star-vs-plus-not-matched` | `heldout/setdiff/GRanges/59` | PB | PB |
| S047 | `e1/granges_core/gr-concat-na-seqlengths-filled-from-right` | `heldout/c/GRanges/57` | PB | PB |
| S048 | `e1/granges_core/seqinfo-merge-conflict-error` | `heldout/merge/Seqinfo/3` | PB | PB |
| S049 | `e1/grid_coverage_na/covna1-partition-num` | `e1/grid_coverage_na/covna9-partition-num-shift` | PB | PB |
| S050 | `heldout/duplicated/GRanges/34` | `heldout/unique/IRanges/8` | AF | PB |
| S051 | `e1/granges_ops/gr-coverage-shift-per-range` | `heldout/coverage/IRanges/41` | AF/PB | PB |
| S052 | `e1/granges_inter/gr-union-conflicting-seqlengths-error` | `heldout/union/GRanges/43` | PB | PB |
| S053 | `e1/grid_mcols_types/mcols-gr-concat-chr-dbl-long` | `e1/grid_mcols_types/mcols-ir-concat-na-lgl-na-chr` | PB | PB |
| S054 | `e1/grid_circular/circ-equal0-count` | `heldout/findOverlaps/GRanges/83` | PB | PB |
| S055 | `e1/grid_iranges_random/rnd0-follow-all` | `heldout/precede/GRanges/50` | PB | PB |
| S056 | `e1/grid_coverage_na/covna1-partition-num-gr` | `e1/grid_coverage_na/covna9-partition-num-gr2` | AF | PB |
| S057 | `e1/granges_parallel/gr-pintersect-zero-width-x` | `heldout/union/IRanges/24` | PB | PB |
| S058 | `e1/grid_bad_args/bad-gr-coverage-shift-dbl` | `heldout/coverage/IRanges/47` | PB | PB |
| S059 | `e1/granges_ops/gr-distance` | `heldout/follow/GRanges/35` | PB | PB |
| S060 | `e1/grid_iranges_random/rnd2-intersect` | `heldout/setdiff/IRanges/4` | PB | PB |
| S061 | `e1/granges_ops/gr-coverage-circular-ignored` | `heldout/coverage/GRanges/97` | PB | PB |
| S062 | `e1/grid_large/large-gr-nested-unique` | `heldout/unique/IRanges/5` | PB | PB |
| S063 | `e1/grid_seqinfo/si0-range` | `heldout/range/GRanges/47` | PB | PB |
| S064 | `e1/construct/ir-concat-mcols-one-sided` | `heldout/c/IRanges/29` | PB | PB |
| S065 | `e1/grid_bad_args/bad-gr-count-maxgap-dbl` | `heldout/subsetByOverlaps/IRanges/16` | PB | PB |
| S066 | `e1/granges_ops/gr-distance-recycle-length-1` | `heldout/distance/IRanges/24` | PB | PB |
| S067 | `e1/grid_iranges_random/rnd0-order` | `heldout/order/IRanges/7` | PB | PB |
| S068 | `e1/granges_ops/gr-fo-seqlength-conflict-error` | `heldout/findOverlaps/GRanges/89` | PB | PB |
| S069 | `e1/grid_bad_args/bad-gr-coverage-shift-NA` | `heldout/gaps/GRanges/51` | PB | PB |
| S070 | `e1/granges_ops/gr-count` | `heldout/countOverlaps/IRanges/3` | PB | PB |
| S071 | `e1/grid_seqinfo/si12-pintersect` | `heldout/pintersect/GRanges/59` | PB | PB |
| S072 | `e1/granges_core/gr-ranges-recycled` | `heldout/GRanges/GRanges/69` | PB | PB |
| S073 | `e1/grid_iranges_random/rnd1-fo-any` | `heldout/findOverlaps/IRanges/58` | PB | PB |
| S074 | `e1/granges_core/gr-mcols-accessor` | `heldout/ranges/GRanges/49` | PB | PB |
| S075 | `e1/grid_iranges_random/rnd0-duplicated` | `heldout/unique/IRanges/25` | PB | PB |
| S076 | `e1/redteam_edges/rt-edge-gr-negative-nearest` | `heldout/precede/GRanges/40` | RA | RA |
| S077 | `e1/granges_core/gr-concat-circular-conflict-error` | `heldout/c/GRanges/38` | PB | PB |
| S078 | `e1/granges_inter/gr-disjoin-empty` | `heldout/reduce/GRanges/94` | PB | PB |
| S079 | `e1/granges_parallel/gr-pintersect-adjacent-pair-is-hit` | `heldout/pintersect/GRanges/54` | PB | PB |
| S080 | `e1/granges_intra/gr-narrow-error` | `heldout/narrow/GRanges/31` | PB | PB |
| S081 | `e1/grid_circular/circ-equal2-nearest-all` | `e1/redteam_circular/rt-circ-nearest-origin-adjacent` | PB | PB |
| S082 | `e1/construct/ir-names-with-na` | `heldout/GRanges/GRanges/67` | PB | PB |
| S083 | `e1/granges_ops/gr-duplicated` | `heldout/unique/GRanges/38` | PB | PB |
| S084 | `e1/grid_construct/ctor-gr-width-lgl` | `e1/grid_construct/ctor-ir-start-width0-width-lgl` | PB | PB |
| S085 | `e1/grid_bad_args/bad-gr-follow-select-NA` | `e1/grid_bad_args/bad-gr-precede-select-neg` | PB | PB |
| S086 | `e1/granges_ops/gr-coverage` | `heldout/coverage/GRanges/96` | PB | PB |
| S087 | `e1/construct/ir-names-na-then-set-mcols` | `heldout/set_names/IRanges/8` | PB | PB |
| S088 | `e1/granges_ops/gr-nearest-empty-subject` | `heldout/precede/GRanges/31` | RA | RA |
| S089 | `e1/grid_bad_args/bad-gr-count-maxgap-lgl` | `e1/grid_bad_args/bad-ir-reduce-mingap-lgl` | PB | PB |
| S090 | `e1/grid_coverage/cov0-gr-plain` | `heldout/coverage/GRanges/84` | PB | PB |
| S091 | `e1/grid_iranges_random/rnd1-union` | `heldout/union/IRanges/8` | PB | PB |
| S092 | `e1/grid_construct/ctor-gr-seq-three-strand-plus` | `heldout/GRanges/GRanges/80` | PB | PB |
| S093 | `e1/granges_intra/gr-promoters-vector` | `heldout/promoters/GRanges/50` | AF | PB |
| S094 | `e1/grid_iranges_random/rnd1-fo-within` | `heldout/findOverlaps/IRanges/13` | PB | PB |
| S095 | `e1/grid_granges/grid-gr-fo-qlong-maxgap` | `heldout/findOverlaps/IRanges/24` | PB | PB |
| S096 | `e1/grid_seqinfo/si25-intersect` | `heldout/intersect/GRanges/55` | PB | PB |
| S097 | `e1/granges_ops/gr-coverage-weight-mcols-column` | `e4/granges-coverage/gr-cov-weight-mcols-name` | AF | AF |
| S098 | `e1/intra/narrow-all-three-error` | `heldout/narrow/IRanges/27` | PB | PB |
| S099 | `e1/redteam_metadata/rt-meta-gr-construct-mcols-divisor` | `e1/redteam_metadata/rt-meta-gr-reserved-width` | PB | PB |
| S100 | `e2/iranges-overlaps/ir.ov.maxgap0.within` | `heldout/findOverlaps/IRanges/27` | PB | PB |
| S101 | `e1/granges_parallel/gr-pintersect-star-x-takes-y-strand` | `heldout/pintersect/GRanges/51` | PB | PB |
| S102 | `heldout/GRanges/GRanges/52` | `heldout/GRanges/GRanges/9` | PB | PB |
| S103 | `e2/granges-core/gr.ctor.factor-seqnames-levels` | `heldout/GRanges/GRanges/50` | PB | PB |
| S104 | `e1/granges_ops/gr-coverage-width-per-seqlevel` | `e4/granges-coverage/gr-cov-width-vector` | AF | PB |
| S105 | `e1/grid_seqinfo/si10-range` | `e1/grid_seqinfo/si37-range` | PB | PB |
| S106 | `e1/granges_ops/gr-follow-empty-range-at-1-second-seqlevel` | `heldout/precede/GRanges/53` | RA | RA |
| S107 | `e1/grid_seqinfo/si10-nearest` | `e1/grid_seqinfo/si37-nearest` | PB | PB |
| S108 | `e1/granges_ops/gr-duplicated-from-last` | `heldout/unique/IRanges/9` | PB | PB |
| S109 | `e1/granges_core/gr-trim-mixed` | `heldout/trim/GRanges/8` | PB | PB |
| S110 | `heldout/intersect/GRanges/43` | `heldout/pintersect/GRanges/58` | PB | PB |
| S111 | `e2/granges-coverage/gr.cov.out-of-bound` | `heldout/coverage/GRanges/92` | PB | PB |
| S112 | `heldout/resize/GRanges/35` | — n=1 | RA | RA |
| S113 | `e1/grid_bad_args/bad-gr-coverage-width-dbl` | `e1/grid_bad_args/bad-ir-coverage-width-dbl` | PB | PB |
| S114 | `e1/grid_bad_args/bad-gr-promoters-up-lgl` | `e1/grid_bad_args/bad-ir-promoters-up-lgl` | PB | PB |
| S115 | `e1/granges_ops/gr-fo-any-maxgap-and-minoverlap-error` | — n=1 | PB | PB |
| S116 | `e1/construct/ir-mcols-accessor-na-names` | — n=1 | RA | PB |
| S117 | `e1/grid_construct/ctor-gr-seq-one-strand-na` | `e2/granges-core/gr.ctor.strand.na` | PB | PB |
| S118 | `heldout/intersect/GRanges/32` | `heldout/intersect/GRanges/51` | PB | PB |
| S119 | `heldout/duplicated/GRanges/36` | `heldout/duplicated/GRanges/40` | PB | PB |
| S120 | `e1/construct/ir-int32-max-ok` | `e4/iranges-edge/ir-edge-width-int-max` | PB | PB |
| S121 | `e1/grid_bad_args/bad-gr-promoters-down-lgl` | `e1/grid_bad_args/bad-ir-promoters-down-lgl` | PB | PB |

## Reproduction / validation

| Command / artifact | Result |
|---|---|
| `/tmp/spotcheck/run.py` (7 existing venv interpreters + Rscript, parallel) | All subprocesses exited0. Initial R input-format mistake subsequently corrected using original pool literals. |
| `Rscript eval/oracle/oracle.R /tmp/spotcheck/r.in.jsonl /tmp/spotcheck/r.out.jsonl` | 201 fresh results; 201/201 equal stored expected projections. |
| `<port-project>/.venv/bin/python -B eval/runner/run_port.py --port <port> --resolved /tmp/spotcheck/<port>.in.jsonl --out /tmp/spotcheck/<port>.out.jsonl` | 235/235 outcomes equal stored status/phase/result/error/reason. |
| `Rscript /tmp/spotcheck/r_symmetry.R` and `Rscript /tmp/spotcheck/r_targeted.R` | Coverage/translation/resize/overflow/empty/NA-name controls above reproduced. |
| BiocPy runner with `expanded.biocpy.in.jsonl`; R oracle with `expanded.r.in.jsonl` | 406 members rerun; fresh R3 universes generated for arbitrary select. |
| `Rscript /tmp/spotcheck/expanded.symmetry.R` | 14/66 S041 and 2/27 S056 order-sensitive; S051 forms split34/2. |
| `PYTHONDONTWRITEBYTECODE=1 eval/ports/biocpy/.venv/bin/python -B /tmp/spotcheck/native.py` | Native API signatures, selected frames, and incomparable-container probe confirmed. |
| `git status --porcelain=v1`; `git diff --quiet`; `git diff --cached --quiet`; six frozen-worktree status/HEAD checks | Clean; no staged files. SHAs: e0 18042e1, e0b 3615722, e1 c5fa8df, e2 70f0c0b, e3 e0e4b41, e4 9e2fa35. |

- No tests added; review-only task. Probe scripts/data exist only in `/tmp/spotcheck`; report is the only durable requested artifact.
- Limits: revised sample is not all200 biocpy clusters; first/last members are not random within clusters; agreement rates are not weighted population estimates; scalar-API and zero-column projection boundaries need consistent policy application before score adjustments.
- Recommended next step: split/co-flag S031, S041, S051, S056; correct S073 evidence; resolve the explicit API-form/projection boundary, then rebuild adjusted summaries. Keep raw scores unchanged pending separate adapter repair and rerun.
