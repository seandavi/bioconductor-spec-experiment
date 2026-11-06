# Re-adjudication verification

**Verdict: accept with notes.** Reviewed `583fc7d` against first adjudication `58424b1`. No repository edits or commits.

## Prior objections: all answered

| Spot-check | Verified disposition |
|---|---|
| S014/S029 | Zero-column GRanges mcols projected with owner length/names; sampled repaired cases pass. |
| S031 | Exact split: **73 pass / 148 port_bug / 3 adapter_fault**. Missing native hits remain wrong (`gr-fo-select-first`); selection projection does not manufacture them. |
| S041 | **22 unsupported / 44 port_bug**. R order-sensitivity no longer excuses the IRanges NA crash. |
| S050 | **18 pass / 19 port_bug** after FALSE→empty-set forwarding repair; independent duplicate-detection failures retained. |
| S051 | **34 unsupported / 2 port_bug**; empty IRanges normalization failures retained. |
| S056/S093/S104 | Documented scalar-only forms now unsupported, not repaired passes. |
| S073 | Corrected to unsorted-query Hits, not countOverlaps dtype. |
| S116 | Explicit port_bug; NA name→`"None"` is not excused as an R accident. |

## Thinness and R accidents

- Audited all adapter/base changes. No search, coverage, ordering, recycling or other native algorithm added. Selected native pairs are only indexed/padded/projected; multiple subjects per query are not silently selected. Zero-column owner context and FALSE→empty container are representation/calling repairs.
- All **108 newly passing case-port pairs** accounted for: 73 selection projections, 18 incomparable-container repairs, 17 zero-column projections. Six frozen port worktrees remain clean at registry SHAs. No score gain found from changing what a port computes.
- Re-ran **every retained r_accident outcome (280)** and checked accident-specific principled answers: 127 shifted-reference nearest-family, 64 empty shapes, 11 local-NA coverages, 8 extreme nearest, 6 resize, 12 average ranks, and 52 remaining validation/overflow/name cases. No separate port defect found hidden by r_accident. Shifted R3 universes were respected. Overflow cases used arithmetic/sequence checks, not an overflowing +1000 reference.

## Required reruns

Ranked baseline clusters by original size where any member's outcome/class changed; thus the 167-member cluster with only two changes is included. Two fresh R/port cases each:

| Baseline cluster | Cases | Confirmed |
|---|---|---|
| e4 coverage/nonlocal NA, 226 | `e1/grid_coverage_na/covna0-partition-int-gr2`; `e1/grid_coverage_na/covna50-empty-inside-int` | port_bug; principled local result/r_accident |
| BiocPy selected Hits, 224 | `e1/granges_ops/gr-fo-empty-subject-select-first`; `heldout/findOverlaps/IRanges/9` | both pass after projection; residual semantic control also rerun |
| BiocPy coverage/type, 167 | `heldout/coverage/GRanges/64`; `heldout/coverage/GRanges/98` | both unsupported vector forms |
| BiocPy coverage/NA crash, 164 | `e1/grid_coverage_na/covna0-partition-int-gr2`; `e1/grid_coverage_na/covna0-partition-int-shift` | unsupported GRanges; port_bug IRanges |
| e0 coverage/nonlocal NA, 148 | `e1/grid_coverage_na/covna0-partition-int-gr2`; `heldout/coverage/IRanges/14` | both unsupported NA representation |

Also reran two members each from the next two clusters (e1/e2, 118 each): port_bug/declared_divergence respectively.

## Notes, not blockers

1. **Three BiocPy empty-query selections still misproject:** `adapters/biocpy.py:312` uses `selected or Hits`; empty RVector is falsy. Already explicitly documented/classified adapter_fault, not port-correct or port_bug. Fix separately with an `is not None` check and regression tests.
2. **Denominator caveat remains material:** 144 former passes become unsupported. Supported/adjusted score movement is partly coverage-scope movement, not improved native fidelity. The revised report discloses this; retain the raw-score/unsupported context.

## Evidence

- `/tmp/readjudication/run.py`: **312/312** fresh port outcomes and **107/107** original R values reproduce stored records; additional 45 shifted R cases and resize control run. All subprocesses exited 0.
- `/tmp/readjudication/verify.py`: all **67,550** scores recompute identically; 665 clusters exhaustively cover adjudicated outcomes, match verdict rules, and reconcile adjusted summaries; all 280 retained R-accident checks pass.
- `PYTHONDONTWRITEBYTECODE=1 eval/.venv/bin/python -B -m pytest -q -p no:cacheprovider eval/tests`: **89 passed**.
- Eval and six frozen worktrees clean; **no staged files**. No tests added (review-only). Scratch probes/results: `/tmp/readjudication/`.

Recommended next step: accept this adjudication; repair the three documented empty-selection adapter outcomes separately.
