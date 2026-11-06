# biocpy adapter notes

Port: `iranges` 0.7.3 + `genomicranges` 0.9.0. Index base: 0-based row indices, coordinates as given (1-based, inclusive).
Header of `biocpy.py` lists the docs used and every translation; this file is the short version.

## Unsupported ops (`unsupported_ops`, 3 of 57)

| Op | Reason |
|---|---|
| `overlapsAny` | no public equivalent |
| `punion`, `psetdiff` | no public parallel union / set difference |

## Unsupported per class or argument (runtime `Unsupported` / `UnknownArgument`)

- `seqlevels<-` maps to the positional `SeqInfo.set_seqnames` only: an unnamed value with as many names as current levels
  (GRanges: `get_seqinfo` -> `set_seqnames` -> `set_seqinfo`). Named values, `pruning.mode`, and add / drop forms are unsupported.
- IRanges has no `rank`, `pintersect`, `ranges`; `GenomicRanges` has no `rank` ties argument (`ties.method`, `na.last`).
- `IRanges(end=)`: the constructor takes `start` / `width` only.
- `GRanges(seqlengths=)`; `GRanges(...)` / `Seqinfo(...)` lacking `seqnames` / `ranges` (the port requires them). The no-argument forms use
  `GenomicRanges.empty()` / `SeqInfo.empty()`.
- `nearest` / `precede` / `follow` with a missing `subject` (R's self-search excludes self; the port has no such form).
- `c(ignore.mcols=, recursive=)`, `range(na.rm=)`, `reduce(with.inframe.attrib=)`, `subsetByOverlaps(invert=)`,
  `findOverlaps(drop.self=, drop.redundant=)`, `seqinfo<-(new2old=, pruning.mode=)`, `use.names` anywhere,
  `ignore.strand` on IRanges / ops whose port method lacks it, `coverage(method=)` on GRanges.
- `[` subscripts whose R and Python meanings differ: zero / negative integers, NA, doubles, logical not of full length.
  Positive integers (shifted to 0-based), full-length logicals and names are passed.
- `order()` of more than one object.
- `findOverlaps` / `countOverlaps` / `subsetByOverlaps` with `type = "equal"` (the port documents any / start / end / within).
- `coverage(weight=, shift=)` given a character column name; `GRanges(seqinfo=<character>)`; `c(<ranges>, <atomic>)`.

## Convention translations

- `duplicated` / `unique` call `biocutils.duplicated` / `biocutils.unique` (the port's public generics; `fromLast` -> `from_last`); results are projected as returned.
- A nested call's ndarray / masked-array result is passed on as an R vector (NA kept), so `x[order(x)]` and `IRanges(start(x), width(x))` work.
- Constructor `...` columns become typed numpy columns (masked when NA is present).

- R `findOverlaps(query, subject)` -> `subject.find_overlaps(query)` (the port's `self` is the subject); same for
  `countOverlaps`, `nearest`, `precede`, `follow`. `subsetByOverlaps(x, ranges)` -> `x.subset_by_overlaps(ranges)`.
- Hits: BiocFrame `query_hits` / `self_hits` (0-based) -> Hits from / to (+1). With `select = "first" | "last" | "arbitrary"`
  the port returns the same frame (one row per matched query); it is projected as R's integer vector of length n_query
  (subject + 1, NA for a query without a row). The rows are not recomputed; a frame that is not one integer pair per query stays Hits.
- 0-based -> 1-based: `order`, nearest / precede / follow integer results (no hit -> NA), `revmap` (only with `with.revmap = TRUE`).
- `[`: positive integer subscripts -> 0-based.
- `seqlengths` / `isCircular` / `genome` accessors: names = seqlevels attached (R returns named vectors).
- Length-1 R vectors -> Python scalars; longer -> `list` (ndarray for `IRanges.coverage` shift / weight);
  named vectors -> dict for the per-seqlevel setters and `gaps(end=)`.
- Missing `subject` in `findOverlaps` / `countOverlaps`: the object is passed as its own query.
- `range(x, y, ...)`: port `combine` first, then `range`. `seqlengths<-` etc. on a GRanges: `get_seqinfo` -> SeqInfo setter -> `set_seqinfo`.
- Coverage: IRanges ndarray -> Rle; GenomicRanges dict -> RleList (names = dict keys, port order).

## Round-trip limits (port, not adapter)

- NA in `names` becomes the string `"None"` in the port (4 distinct inputs).
- Ranges whose end exceeds int32 are rejected at construction (build-phase port error, 6 inputs).
