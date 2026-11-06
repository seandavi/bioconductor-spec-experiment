# e0b adapter notes

Port: `ranges_e0b` (frozen). Index base: 0-based positional indices, `NO_HIT = -1` for NA; coordinates 1-based.

## Unsupported ops (4 of 57 in ops.json)

`seqlevels<-`, `seqlengths<-`, `isCircular<-`, `genome<-`: the port has no setters, only `GRanges.with_seqinfo`.
The other 53 ops have an `op_*`.

## Unsupported arguments / forms (per case, `unsupported`)

- Arguments the port lacks: `use.names` (all ops), `with.revmap` on range/disjoin, `fromLast`, `incomparables`,
  `na.last`, `na.rm`, `ignore.mcols`, `recursive`, `drop.self`, `drop.redundant`, `drop.nohit.ranges`,
  `strict.strand`, `new2old`, `pruning.mode`, `with.inframe.attrib`, `nthread`, `method`, `by`, GRanges `seqlengths`, `j`/`drop` on `[`.
- `range(x, y, ...)`, `merge(x, y, z)`, multi-key `order(a, b)`: the port takes one object.
- Subscripts: negative/zero/NA integers (port: negative counts from the end), double, NA logical,
  logical recycling (port: mask of exactly length n), NULL.
- NA inside a numeric/logical vector argument (the port has no NA in int/bool arrays).
- A factor as GRanges `seqnames` (the port has no factor type).

## Convention translations

- Names `a.b` -> `a_b`; `y` -> `other`, `ranges` (subsetByOverlaps) -> `subject`; generics -> methods.
- 1-length vector -> scalar, longer -> ndarray / list; named vector `start`/`end`/`width` -> `Mapping[str, int]` (GRanges.gaps/coverage, documented).
- Subscript: positive integers `- 1`. Results `+ 1`: Hits, `select=` results, `order`, `reduce` revmap. `NO_HIT` -> NA.
- `distance`: `-1` -> NA (documented as R NA in `GRanges.distance`).
- mcols: dict of polars Series (int/dbl/chr/lgl -> Int64/Float64/String/Boolean, Factor -> Enum, IntegerList -> List(Int64)).
- `reduce(with_revmap=True)` tuple -> result with a `revmap` IntegerList mcols column.
- `seqlengths`/`isCircular`/`genome`: GRanges dicts and Seqinfo tuples -> vectors named by seqlevel.
- GRanges coverage `dict[str, Rle]` -> RleList; `mcols` zero-column frame projects with its own height (0).
