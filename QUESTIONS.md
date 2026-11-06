# Open questions for Bioconductor developers

Six agent-built ports were told that R is normative. Where R's behavior
looks like an implementation accident rather than a contract, "R is
normative" makes every port copy the accident or record a divergence. These
are the cases where a person who knows the code should decide. Each has a
reproducible R snippet (R 4.6.0, IRanges 2.47.2, GenomicRanges 1.65.1).

## 1. `coverage()` with NA weights: is NA local, or does it propagate?

```r
x <- IRanges(c(1, 11, 20), width = c(5, 3, 9))
w <- c(1L, NA, 2L)
as.vector(coverage(x, weight = w))                    # auto (sort/hash here)
#  1 1 1 1 1 0 0 0 0 0 NA NA NA NA NA NA NA NA NA ... NA   (positions 11-28)
as.vector(coverage(x, weight = w, method = "naive"))
#  1 1 1 1 1 0 0 0 0 0 NA NA NA 0 0 0 0 0 0 2 2 2 2 2 2 2 2 2
```

The cumulative-sum methods carry the NA to every later position, including
positions only the weight-2 range covers; `naive` keeps it local. Which one
is the contract? 1,239 of the 9,650 evaluation cases turn on this. The ports
split: E1 reverse-engineered and cloned the propagating behavior, E2 declared
a divergence, E4 and BiocPy count NA as 0.

## 2. Seqlevel order when no seqinfo is given

```r
seqlevels(GRanges(c("chr2", "chr10", "chr1", "chrX"), IRanges(1, width = 1)))
#  "chr2" "chr10" "chr1" "chrX"     (first appearance)
```

BiocPy genomicranges 0.9.0 sorts lexically (`chr1 chr10 chr2 chrX`); all six
ports follow R. Is first-appearance order normative, so that other
implementations should change, or implementation-defined, so that
comparisons should ignore order unless seqinfo is explicit?

## 3. `rank()` on IRanges without `ties.method`

```r
rank(IRanges(c(5, 1, 5), width = 1))
#  Error: 'arg' must be of length 1
```

This looks like a `match.arg()` slip rather than a design choice. Every port
defaults to `"average"` like `base::rank()` and records the divergence.
Should R's default work?

## 4. Must ports reject what R rejects?

```r
reduce(GRanges("chr1", IRanges(1:2, width = 3)), ignore.strand = NA)
#  Error: 'ignore.strand' must be TRUE or FALSE
```

Every port except E1 silently accepts invalid flags like this. It's the most
common defect shared across ports, and it's invisible unless the comparison
includes inputs R refuses. Is argument validation part of the contract a
port should match, or only the results on valid input?

## 5. Edge behavior near int32 limits and coordinates ≤ 0

E4 (ADR-0010) found `gaps`, `setdiff` and `disjoin` misbehaving at
`end = 2^31 - 1`, and a non-monotonic `maxgap` threshold near the limit (a
hit at `maxgap = 2^31 - 3` disappears at `2^31 - 2`). On GRanges,
`precede`/`follow`/`nearest` with start or end ≤ 0 can leak across
sequences, apparently because sequences share one number line with
sentinels. Are these worth fixing upstream, and should ports reproduce them
in the meantime? The ports declined.

## 6. Is circular-sequence semantics in scope for a port?

`isCircular()` changes overlap, nearest and coverage results. Only E1 (via
its red team) implemented it; E2 and E4 silently compute linear answers;
E0, E0b and E3 declared it out of scope. Is a port that ignores circularity
acceptable if it says so?

Evidence for every item is in `eval/adjudication/` (clusters and probes) and
`ledger/LEDGER.md`.
