# Comparator

R is normative. Comparison is **exact by default**. This file is the only
place a relaxation may be defined; there are exactly three.

## Outcomes per (case, port)

Evaluated top to bottom; first match wins.

| # | R (expected) | Port record | Outcome |
|---|---|---|---|
| 0 | `excluded` | — | not scored (listed in report with reason) |
| 1 | any | `unsupported` | `unsupported` |
| 2 | any | `adapter_error` | `adapter_error` |
| 2b | any | `port_error`, phase `build` | `port_error_unexpected` (R built that input; a port that cannot gets no R1 pass) |
| 2c | any | `port_error`, phase `project` | `port_error_unexpected` (the operation returned; failing to read its result back is not an error from the operation, so no R1 pass) |
| 3 | error | `port_error`, phase `call`, not a timeout | `pass` (relaxation R1) |
| 4 | error | value | `mismatch` (diff `r_error_port_value`) |
| 4b | error | `port_error`, origin `timeout` | `port_error_unexpected` (a hang is not an error from the operation) |
| 5 | value | `port_error` (incl. timeout) | `port_error_unexpected` |
| 6 | value | value | compare (rule below) → `pass` \| `mismatch` |

Port record statuses are assigned by the runner (ADAPTERS.md, "Classification").
`unsupported` and `adapter_error` never earn a pass, even when R errored.

## Value rules (row 6)

| Rule | Applies when | Pass iff |
|---|---|---|
| R3 `arbitrary_in_universe` | `findOverlaps` with `kwargs.select == "arbitrary"`, or `nearest` with `select` absent or `"arbitrary"`; and R result kind is `vector` | port is `vector` of R's `type`, same length, `names` equal; per query *i*: allowed = `universe` hits of *i* ∪ {R's own value for *i*, if not NA}; port[*i*] is `null` iff allowed is empty, else port[*i*] ∈ allowed |
| R2 `hits_query_set` | R result is `Hits` and `sorted_by_query` is true | port is `Hits`; `nLnode`, `nRnode` equal; port `from` non-decreasing; sorted (from, to) pair lists equal (multiset) |
| `exact` | everything else | deep equality (below) |

Exact equality: same `kind`; same keys; lists equal in length and order;
int/chr/lgl leaves equal; dbl leaves equal as IEEE doubles (`"NaN"` = `"NaN"`,
`"Inf"` = `"Inf"`, `-0` = `0`); `null` equals only `null`.

## Relaxations and why

| Id | Relaxation | Justification |
|---|---|---|
| R1 | both error → pass; messages not compared | Error text is R-specific; error *occurrence* is the contract. |
| R2 | Hits compared as a set within each query (query order still exact) | R returns `SortedByQueryHits`: sorted query order is a class invariant, subject order within a query is NCList traversal order and undocumented. Measured in R 4.6 / IRanges 2.47.2: 403 of 3000 random GRanges `findOverlaps` calls emitted within-query subject indices out of order. Applies to every R result whose class is `SortedByQuery*` (`findOverlaps`, `nearest/precede/follow(select="all")`), decided by the oracle's `sorted_by_query`, not by op name. |
| R3 | `select = "arbitrary"`: any of R's `select = "all"` hits | `?findOverlaps`: "arbitrary" returns an arbitrary overlapping subject. `?nearest`: ties are resolved "arbitrary for nearest" (its default). The oracle stores the `select = "all"` result of the same call as `universe`. R's own pick is always allowed too: in R 4.6 / GenomicRanges 1.65.1, `nearest(<GRanges>)` picks the lower-index tie while `select="all"` returns only the higher one (`nearest(GRanges("chr1",IRanges(20,21)), GRanges("chr1",IRanges(c(10,27,33),width=c(5,3,1))))` is 1, `select="all"` is {2}); IRanges returns both. |

Not relaxed (deliberately): doubles (no tolerance; see diff kind `float_ulp`),
vector `names`, storage `type`, row order of ranges, `revmap` element order,
`precede`/`follow` default `select` (documented first/last).

## Mismatch detail (informational, never changes the outcome)

Each mismatch records `path` (first differing JSON path), `r` and `port`
values there (truncated), and `kinds`: the set over **all** differing leaves:

| Kind | Meaning |
|---|---|
| `kind` | different result kinds |
| `length` | list lengths differ |
| `type` | vector/Rle storage type differs |
| `names` | only `names` / `rownames` differ |
| `value` | a leaf value differs |
| `float_ulp` | every differing leaf is dbl within 4 ulp |
| `order` | same multiset, different order (or R2's `from` not sorted) |
| `r_error_port_value` | row 4 |
| `not_in_universe` | R3 failed |

## Scores

`pass rate = pass / scored` (all five outcomes in the denominator). Reported
per port x source and per port x family x source; heldout is the headline
(no port has seen it). Excluded and unconvertible cases are counted, never dropped silently.
