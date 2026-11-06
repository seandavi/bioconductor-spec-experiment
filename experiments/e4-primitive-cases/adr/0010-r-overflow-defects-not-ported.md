# ADR-0010: R integer-overflow defects are not ported

- Status: accepted. Divergence from R (deliberate); cases carry `python: {skip: ...}`.
- Where the R golden is an artefact of 32-bit overflow or an inconsistent
  return type, Python does the sensible thing and the case is skipped (not
  inverted), so the golden still documents R's output.

| Case | R | Python |
|---|---|---|
| `ir-edge-disjoin-end-int-max` | error computing `end+1` at 2^31-1 | normal result |
| `ir-edge-fo-maxgap-int-max` | `maxgap=2^31-1` overflows, no hit | finds the hit |
| `gr-edge-empty-precede-returns-hits-in-r` | zero-length query returns a `Hits` | select-shaped (empty masked) vector |
| `gr-edge-empty-follow-returns-hits-in-r` | same | same |

- Python behaviour for the first two is covered by `tests/test_edge.py`.

Added after probing (all end at 2^31-1):

| Case | R | Python |
|---|---|---|
| `ir-edge-gaps-to-max` | `gaps(end=2^31-1)` returns `[1,2^31-1],[8,2^31-1]` (wrong) | `[1,4],[8,2^31-1]` |
| `ir-edge-setdiff-max` | `setdiff(1..2^31-1, 5..7)` returns `[1,2^31-1]` (hole lost) | `[1,4],[8,2^31-1]` |
