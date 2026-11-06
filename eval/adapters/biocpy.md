# biocpy adapter notes

Port: `iranges` 0.7.3 + `genomicranges` 0.9.0. Index base: 0-based row indices, coordinates as given (1-based, inclusive).
Header of `biocpy.py` lists the docs used and every translation; this file is the short version.

## Unsupported ops (`unsupported_ops`, 6 of 59)

| Op | Reason |
|---|---|
| `seqlevels<-` | no public setter with R's match-by-name / drop / add semantics (`SeqInfo.set_seqnames` renames positionally) |
| `overlapsAny` | no public equivalent |
| `punion`, `psetdiff` | no public parallel union / set difference |
| `duplicated`, `unique` | no public equivalent |

## Unsupported per class or argument (runtime `Unsupported` / `UnknownArgument`)

- IRanges has no `rank`, `pintersect`, `ranges`; `GenomicRanges` has no `rank` ties argument (`ties.method`, `na.last`).
- `IRanges(end=)`: the constructor takes `start` / `width` only.
- `GRanges(seqlengths=)`, `GRanges()` / `Seqinfo()` with no `seqnames` / `ranges` (the port requires them).
- `nearest` / `precede` / `follow` with a missing `subject` (R's self-search excludes self; the port has no such form).
- `c(ignore.mcols=, recursive=)`, `range(na.rm=)`, `reduce(with.inframe.attrib=)`, `subsetByOverlaps(invert=)`,
  `findOverlaps(drop.self=, drop.redundant=)`, `seqinfo<-(new2old=, pruning.mode=)`, `use.names` anywhere,
  `ignore.strand` on IRanges / ops whose port method lacks it, `coverage(method=)` on GRanges.
- `[` subscripts whose R and Python meanings differ: zero / negative integers, NA, doubles, logical not of full length.
  Positive integers (shifted to 0-based), full-length logicals and names are passed.
- `order()` of more than one object.

## Convention translations

- R `findOverlaps(query, subject)` -> `subject.find_overlaps(query)` (the port's `self` is the subject); same for
  `countOverlaps`, `nearest`, `precede`, `follow`. `subsetByOverlaps(x, ranges)` -> `x.subset_by_overlaps(ranges)`.
- Hits: BiocFrame `query_hits` / `self_hits` (0-based) -> Hits from / to (+1); `select = "first" | "last" | "arbitrary"`
  in findOverlaps also returns that BiocFrame in the port, so it is projected as Hits (not an integer vector).
- 0-based -> 1-based: `order`, nearest / precede / follow integer results (no hit -> NA), `revmap` (only with `with.revmap = TRUE`).
- `[`: positive integer subscripts -> 0-based.
- `seqlengths` / `isCircular` / `genome` accessors: names = seqlevels attached (R returns named vectors).
- Length-1 R vectors -> Python scalars; longer -> `list` (ndarray for `IRanges.coverage` shift / weight);
  named vectors -> dict for the per-seqlevel setters and `gaps(end=)`.
- Missing `subject` in `findOverlaps` / `countOverlaps`: the object is passed as its own query.
- `range(x, y, ...)`: port `combine` first, then `range`. `seqlengths<-` etc. on a GRanges: `get_seqinfo` -> SeqInfo setter -> `set_seqinfo`.
- Coverage: IRanges ndarray -> Rle; GenomicRanges dict -> RleList (names = dict keys, port order).

## Round-trip limits (port, not adapter)

- NA in `names` becomes the string `"None"` in the port (11 inputs).
- Ranges whose end exceeds int32 are rejected at construction (build-phase port error).
