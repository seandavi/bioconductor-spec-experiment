Implemented 119 adversarial cases in four commits. **67 agree; 52 strict xfails** expose 47 bug cases and 5 ADR judgments. No `src/` edits or existing expectation changes.

Complete per-case tables: [circular](experiments/e1-differential/redteam/circular.md), [metadata](experiments/e1-differential/redteam/metadata.md), [edges](experiments/e1-differential/redteam/edges.md).

| Case id | R | Python | View |
|---|---|---|---|
| `rt-circ-origin-any` | Hit across circular origin | No hit | Bug |
| `rt-meta-ir-set-mcols-scalar` | Recycles metadata rows | Errors | Bug |
| `rt-meta-gr-pintersect-existing-hit-keep` | Replaces `hit` column | Duplicate-column error | Bug |
| `rt-edge-ir-default-empty` | Empty IRanges | Errors | Bug |
| `rt-meta-gr-set-seqinfo-reorder` | Requires explicit mapping | Remaps implicitly | ADR judgment |
| `rt-edge-gr-negative-follow` | Errors | Returns geometric hits | ADR judgment |

Changed files: cases, caches, generators, reports, and harness support only.

Validation: **3888 passed, 52 xfailed**. Working tree clean; nothing staged.

Open risks/questions: divergences remain unresolved; independent review required.  
Recommended next step: review findings, then fix Python or approve targeted ADRs.