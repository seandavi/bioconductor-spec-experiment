# Target: IRanges / GRanges in Python

Every experiment builds an **independent Python implementation** of the
range classes below, with R/Bioconductor as the normative reference
(IRanges 2.47.x, GenomicRanges 1.65.x, Bioc 3.24 — whatever `Rscript`
on this machine reports). Do not depend on, vendor, or copy BiocPy
(`iranges`, `genomicranges`, `biocframe` on PyPI); the point is an
independent port.

Where R and your implementation differ, R is right unless a recorded
decision (in whatever form your experiment uses) says otherwise.

## Required scope

### IRanges (integer ranges, 1-based, closed: `end = start + width - 1`)

- Construction from any two of start/end/width; names; zero-width ranges
- Accessors: `start`, `end`, `width`, `names`, `length`
- Subsetting (integer, logical, names), concatenation
- Mcols (per-range metadata columns)
- Intra-range: `shift`, `narrow`, `resize` (fix = start/end/center),
  `flank` (start, both), `promoters`
- Inter-range: `range`, `reduce` (incl. `min.gapwidth`,
  `drop.empty.ranges`, `with.revmap`), `gaps` (incl. `start`/`end`),
  `disjoin`
- Set ops: `union`, `intersect`, `setdiff` (and parallel `p*` versions)
- Overlaps: `findOverlaps` (type = any/start/end/within/equal,
  `maxgap`, `minoverlap`, `select` = all/first/last/arbitrary),
  `countOverlaps`, `overlapsAny`, `subsetByOverlaps`
- Nearest: `nearest`, `precede`, `follow`, `distance`
- `coverage` (returning run-length data in some form)
- Ordering: `order`, `sort`, `rank`, `duplicated`, `unique`

### GRanges (IRanges + seqnames + strand + seqinfo)

- Construction; `seqnames`, `strand` (`+`, `-`, `*`), `ranges`, mcols
- Seqinfo: `seqlevels`, `seqlengths`, `isCircular`, `genome`;
  merging seqinfo on combine; `trim` against seqlengths
- All IRanges operations above, strand-aware where R is strand-aware,
  including `ignore.strand`
- Ordering semantics across seqlevels and strand

## Out of scope

GRangesList / IRangesList, Rle as a general class (internal use is
fine), Views, genome annotation lookups, file I/O (BED/GFF), plotting.

## Done criteria

An experiment is done when:

1. Every required operation is implemented as a public API.
2. The experiment's own test suite passes (`uv run pytest`).
3. The experiment's own README describes its design and how its
   artifacts (if any: specs, ADRs, conformance cases, harnesses) are
   meant to be maintained.
4. All work is committed in stages (see PROTOCOL.md).
