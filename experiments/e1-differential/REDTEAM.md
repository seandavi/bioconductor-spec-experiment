# Red-team results

119 new differential cases: **67 agree, 52 strict expected failures**.
47 failing cases look like bugs; 5 need ADR judgment. No difference was accepted.
`src/` and pre-existing test expectations were not changed.

## Per-case tables

Every row gives the case id, actual R result, actual Python result, and verdict.
Outputs are abbreviated; committed caches retain complete R values.

| Batch | Cases | Agree | Xfail | Table |
|---|---:|---:|---:|---|
| Circular seqinfo | 33 | 16 | 17 | [circular](redteam/circular.md) |
| Metadata / validation | 40 | 15 | 25 | [metadata](redteam/metadata.md) |
| Empty / ordering / boundaries | 46 | 36 | 10 | [edges](redteam/edges.md) |

## Findings

| Behavior | Failing cases | View |
|---|---:|---|
| Circular overlaps omit origin/periodic hits; nearest inherits it | 17 | Bug |
| Short mcols rows do not recycle | 12 | Bug |
| Masked logical flank arguments lose NA | 2 | Bug |
| Coverage cannot resolve some metadata-column arguments | 2 | Bug |
| Coverage rejects valid missing numeric weights | 5 | Bug |
| pintersect appends an existing `hit` column instead of replacing it | 1 | Bug |
| Literal strand `"NA"` becomes `"*"` | 1 | Bug |
| Seqinfo accepts nonlogical circularity / numeric names | 3 | Bug |
| Recycled GRanges coordinates retain unrecycled names | 1 | Bug |
| Default empty constructors error | 2 | Bug |
| Blank-name subscript selects instead of erroring | 1 | Bug |
| Seqinfo setter silently remaps/drops unused levels | 2 | ADR judgment |
| Negative genomic coordinates break R nearest-family sentinel arithmetic | 3 | ADR judgment |

The last group includes an R nearest hit across different sequences and an R follow error.
Do not blindly copy those results. Decide whether to mirror them or record a divergence.
The seqinfo setter may reasonably be more convenient in Python, but needs an explicit decision.

## Validation / maintenance

- R 4.6.0; IRanges 2.47.2; GenomicRanges 1.65.1; S4Vectors 0.51.6.
- Baseline: 3808 passed. Final full suite: **3888 passed, 52 xfailed**.
- Fresh R run of new cases: **67 passed, 52 xfailed**.
- Added 13 harness regression tests. Fixed one harness-only issue: typed empty
  numpy metadata lost its dtype through `tolist()`. Four apparent port failures disappeared.
- Xfail has a reason, is strict, and catches comparison failures only. Harness crashes still fail.
- Generators preserve reviewed xfail reasons. Reports refuse stale xfail status.

```sh
uv run pytest tests/differential -k 'rt-circ or rt-meta or rt-edge' --refresh-r
RANGES_E1_NO_R=1 uv run pytest
uv run python tests/differential/gen/redteam_report.py
# Inspect any unresolved case without xfail masking:
uv run pytest tests/differential -k rt-circ-origin-any --runxfail
```

Judgment calls and the corrected timestamp are in append-only `NOTES.md`.
Independent review remains required. Fix Python or approve a short ADR before removing xfail marks.
