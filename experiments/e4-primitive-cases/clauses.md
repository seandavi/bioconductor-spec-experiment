# Clauses

One behavior per row. Append-only; retire with `(retired)`.
Every clause needs >=1 case (checked by `tests/test_case_inventory.py`).

| id | statement | R ref |
|---|---|---|
| IR-CON-START-WIDTH-1 | `IRanges(start, width)`: `end = start + width - 1` | `?IRanges` |
| IR-CON-START-END-1 | `IRanges(start, end)`: `width = end - start + 1` | `?IRanges` |
| IR-CON-END-WIDTH-1 | `IRanges(end, width)`: `start = end - width + 1` | `?IRanges` |
| IR-CON-EMPTY-1 | zero-length IRanges is valid; all accessors return empty vectors | `?IRanges` |
| IR-CON-ZEROWIDTH-1 | zero-width ranges are valid: `end = start - 1` | `?IRanges` |
| IR-CON-NEGWIDTH-1 | negative width is an error | `?IRanges` |
| IR-CON-NAMES-1 | names are carried through; NA names allowed | `?IRanges` |
| IR-CON-MCOLS-1 | mcols columns are carried through, in order | `?mcols` |
| IR-CON-RECYCLE-1 | length-1 `start`/`end`/`width` recycle to the longest | `?IRanges` |
| IR-CON-ALL3-1 | all three of start/end/width allowed when consistent; inconsistent is an error | `?IRanges` |
| IR-CON-TWO-1 | fewer than two of start/end/width is an error | `?IRanges` |
| IR-CON-ENDLT-1 | `end < start - 1` is an error | `?IRanges` |
| IR-CON-INT32-1 | coordinates beyond int32: R yields NA end with a warning; Python raises (ADR-0006) | `?IRanges` |
| IR-ACC-1 | `start`, `end`, `width`, `names`, `length` return the stored values | `?IRanges` |
| IR-SUB-INT-1 | integer subscript selects in the given order, repeats allowed | `?"[,Vector-method"` |
| IR-SUB-LGL-1 | logical subscript selects TRUE positions | `?"[,Vector-method"` |
| IR-SUB-NAMES-1 | character subscript selects by name (first match) | `?"[,Vector-method"` |
| IR-SUB-ERR-1 | out-of-bounds index or unknown name is an error | `?"[,Vector-method"` |
| IR-SUB-KEEP-1 | subsetting keeps names and mcols aligned | `?"[,Vector-method"` |
| IR-CAT-1 | `c()` concatenates in argument order | `?"c,Vector-method"` |
| IR-CAT-NAMES-1 | `c()` pads unnamed pieces with `""` when any piece has names | `?"c,Vector-method"` |
| IR-CAT-MCOLS-1 | `c()` unions mcols columns; missing cells are NA | `?"c,Vector-method"` |
| IR-CAT-EMPTY-1 | `c()` with a zero-length piece is a no-op for that piece | `?"c,Vector-method"` |
| IR-SHIFT-1 | `shift` adds the offset to start; width unchanged; names/mcols kept | `?"intra-range-methods"` |
| IR-SHIFT-2 | `shift` accepts a per-range vector | `?"intra-range-methods"` |
| IR-SHIFT-3 | `shift` with a vector that recycles cyclically: R accepts, Python errors (ADR-0007) | `?"intra-range-methods"` |
| IR-NARROW-1 | `narrow(start=)` / `narrow(end=)`: positive counts from range start, negative from range end | `?narrow` |
| IR-NARROW-2 | `narrow` with start+width, end+width, or start+end | `?narrow` |
| IR-NARROW-3 | `narrow` may produce an empty range at `start = width+1` | `?narrow` |
| IR-NARROW-4 | `narrow` that would widen, give negative width, or supply 3 (or only width) args is an error | `?narrow` |
| IR-RESIZE-1 | `resize(fix="start")` keeps start; `"end"` keeps end | `?resize` |
| IR-RESIZE-2 | `resize(fix="center")` shifts start by `floor((oldwidth - width)/2)` | `?resize` |
| IR-RESIZE-3 | `resize` accepts per-range width and fix; `use.names=FALSE` drops names | `?resize` |
| IR-RESIZE-4 | `resize` to width 0 is valid; negative width or unknown `fix` is an error | `?resize` |
| IR-FLANK-1 | `flank(width>0)` is the region just before start (`start=TRUE`) or after end | `?flank` |
| IR-FLANK-2 | `flank(width<0)` is the region inside the range at the chosen side | `?flank` |
| IR-FLANK-3 | `flank(both=TRUE)` spans `abs(width)` each side of the boundary | `?flank` |
| IR-FLANK-4 | `flank(width=0)` gives an empty range at the boundary | `?flank` |
| IR-PROMOTERS-1 | `promoters`: `[start-upstream, start+downstream-1]` (defaults 2000/200) | `?promoters` |
| IR-PROMOTERS-2 | `promoters` with negative upstream/downstream is an error | `?promoters` |
| IR-INTRA-EMPTY-1 | intra-range ops on a zero-length IRanges return zero-length IRanges | `?"intra-range-methods"` |
| IR-RANGE-1 | `range` is `[min(start), max(end)]` over all ranges (zero-width included); drops names/mcols | `?range` |
| IR-RANGE-2 | `range` of a zero-length IRanges is zero-length; `range(x, y)` pools both | `?range` |
| IR-REDUCE-1 | `reduce` merges overlapping ranges; adjacent ranges merge at the default `min.gapwidth=1` | `?reduce` |
| IR-REDUCE-2 | `min.gapwidth=0` keeps adjacent ranges separate; larger values also merge gaps `< min.gapwidth` | `?reduce` |
| IR-REDUCE-3 | zero-width ranges are kept unless they touch/lie inside another range (default) | `?reduce` |
| IR-REDUCE-4 | `drop.empty.ranges=TRUE` removes zero-width ranges before merging | `?reduce` |
| IR-REDUCE-5 | `reduce` result has no names or input mcols; input order and duplicates do not matter | `?reduce` |
| IR-REDUCE-6 | `with.revmap=TRUE` adds `revmap`: member indices per result range, in `order(x)` order | `?reduce` |
| IR-REDUCE-7 | negative `min.gapwidth` is an error; zero-length input reduces to zero-length | `?reduce` |
| IR-GAPS-1 | `gaps` default bounds are the extent of the non-empty coverage: only interior gaps | `?gaps` |
| IR-GAPS-2 | `gaps(start=, end=)` complements within `[start, end]`, including leading/trailing gaps | `?gaps` |
| IR-GAPS-3 | `gaps` with one bound given uses the coverage extent for the other | `?gaps` |
| IR-GAPS-4 | `gaps` ignores zero-width ranges; no non-empty coverage and a missing bound gives zero-length | `?gaps` |
| IR-GAPS-5 | `gaps` with `end < start - 1` or all-covered bounds gives zero-length | `?gaps` |
| IR-DISJOIN-1 | `disjoin` splits at every start and end+1; adjacent ranges stay separate | `?disjoin` |
| IR-DISJOIN-2 | `disjoin` ignores zero-width ranges | `?disjoin` |
| IR-DISJOIN-3 | `disjoin(with.revmap=TRUE)` lists covering range indices per segment, ascending | `?disjoin` |
| IR-DISJOIN-4 | `disjoin` of zero-length input is zero-length; duplicates collapse | `?disjoin` |
| IR-UNION-1 | `union` = non-empty coverage of x and y, normalised (sorted, merged incl. adjacent) | `?"setops"` |
| IR-UNION-2 | `union` with zero-length input | `?"setops"` |
| IR-UNION-3 | `union` merges adjacent ranges | `?"setops"` |
| IR-UNION-4 | `union` drops zero-width ranges | `?"setops"` |
| IR-UNION-5 | `union` normalises x even when y is empty; result has no names/mcols | `?"setops"` |
| IR-INTERSECT-1 | `intersect` = positions covered by both | `?"setops"` |
| IR-INTERSECT-2 | `intersect` with zero-length input is zero-length | `?"setops"` |
| IR-INTERSECT-3 | `intersect` of adjacent ranges is empty | `?"setops"` |
| IR-INTERSECT-4 | `intersect` drops zero-width ranges | `?"setops"` |
| IR-INTERSECT-5 | `intersect` works on coverage, not on individual ranges | `?"setops"` |
| IR-SETDIFF-1 | `setdiff(x, y)` = covered by x and not by y | `?"setops"` |
| IR-SETDIFF-2 | `setdiff` with zero-length input | `?"setops"` |
| IR-SETDIFF-3 | `setdiff` may split or remove ranges | `?"setops"` |
| IR-SETDIFF-4 | `setdiff` ignores zero-width ranges in y; drops them in x | `?"setops"` |
| IR-PUNION-1 | `punion` = `[min start, max end]` per pair; a gap between the pair is an error | `?"setops"` |
| IR-PUNION-2 | `punion(fill.gap=TRUE)` fills the gap | `?"setops"` |
| IR-PUNION-3 | `punion` gap test uses `max(start) > min(end) + 1`, so a zero-width range beside the other member is fine | `?"setops"` |
| IR-PUNION-4 | parallel ops require equal lengths (no recycling) | `?"setops"` |
| IR-PINTERSECT-1 | `pintersect` of a pair with a gap is an ambiguous empty range: error by default | `?"setops"` |
| IR-PINTERSECT-2 | `pintersect` result keeps x's names, drops mcols | `?"setops"` |
| IR-PINTERSECT-3 | adjacent or zero-width pairs give an unambiguous empty range at `[max start, min end]` | `?"setops"` |
| IR-PINTERSECT-4 | `resolve.empty="start"` puts ambiguous empties at start(x); `"max.start"` at max start | `?"setops"` |
| IR-PSETDIFF-1 | `psetdiff` subtracts y from x pairwise | `?"setops"` |
| IR-PSETDIFF-2 | `psetdiff` outcomes by overlap geometry (left, right, cover, disjoint) | `?"setops"` |
| IR-PSETDIFF-3 | zero-width x is returned unchanged | `?"setops"` |
| IR-PSETDIFF-4 | y strictly inside x (would split it) is an error, including zero-width y | `?"setops"` |
| IR-PSETDIFF-5 | zero-width y at x's start edge or beyond its end is not an error | `?"setops"` |
| IR-SETDIFF-5 | `setdiff` handles nested/multiple overlaps on coverage | `?"setops"` |
| IR-FO-ANY-1 | `findOverlaps` default (`type="any"`, `maxgap=-1`) hits ranges sharing at least one position | `?findOverlaps` |
| IR-FO-ANY-2 | adjacent ranges (gap 0) do not overlap by default | `?findOverlaps` |
| IR-FO-ANY-3 | zero-width range hits only a range strictly containing its position (needs a position on both sides); `maxgap=0` also hits boundaries; `minoverlap>=1` never hits | `?findOverlaps` |
| IR-FO-ANY-4 | `maxgap=g`: hit iff gap `<= g` (gap = distance between nearest ends, 0 when adjacent) | `?findOverlaps` |
| IR-FO-ANY-5 | `minoverlap=m`: hit iff shared width `>= m` | `?findOverlaps` |
| IR-FO-ERR-1 | `maxgap < -1`, `minoverlap < 0`, both non-default with type any, bad `type`/`select` are errors | `?findOverlaps` |
| IR-FO-START-1 | `type="start"` (default maxgap 0): starts equal | `?findOverlaps` |
| IR-FO-START-2 | `type="start"`, `maxgap=g`: `abs(start diff) <= g`; ranges need not overlap | `?findOverlaps` |
| IR-FO-END-1 | `type="end"`: ends equal | `?findOverlaps` |
| IR-FO-END-2 | `type="end"`, `maxgap=g`: `abs(end diff) <= g` | `?findOverlaps` |
| IR-FO-END-3 | start/end types with zero-width queries compare coordinates (`[s, s-1]`) | `?findOverlaps` |
| IR-FO-WITHIN-1 | `type="within"`: query inside subject (`ss <= qs`, `qe <= se`) | `?findOverlaps` |
| IR-FO-WITHIN-2 | `within` + `maxgap=g>0`: also `subject width - query width <= g`; `maxgap` 0 or -1 means no limit (R quirk) | `?findOverlaps` |
| IR-FO-WITHIN-3 | zero-width query within a subject includes its boundaries; zero-width subject contains no non-empty query | `?findOverlaps` |
| IR-FO-EQUAL-1 | `type="equal"` (default maxgap 0): same start and end | `?findOverlaps` |
| IR-FO-EQUAL-2 | `equal` + `maxgap=g`: both start and end differences `<= g` | `?findOverlaps` |
| IR-FO-EQUAL-3 | zero-width equal compares `[s, s-1]` | `?findOverlaps` |
| IR-FO-MINOV-1 | `minoverlap` also applies to start/end/within/equal types | `?findOverlaps` |
| IR-FO-ORDER-1 | all-hits are sorted by query then subject; Hits carry query/subject lengths | `?findOverlaps` |
| IR-FO-ORDER-2 | `findOverlaps(x)` with no subject uses x as subject (self hits included) | `?findOverlaps` |
| IR-FO-EMPTY-1 | zero-length query or subject gives zero hits; lengths preserved | `?findOverlaps` |
| IR-FO-SELECT-1 | `select="first"`: lowest subject index per query; NA when no hit | `?findOverlaps` |
| IR-FO-SELECT-2 | `select="last"`: highest subject index per query; NA when no hit | `?findOverlaps` |
| IR-FO-SELECT-3 | `select="arbitrary"`: some hit per query (implementation-defined; ADR-0005) | `?findOverlaps` |
| IR-FO-SELECT-4 | `select` with zero-length query or subject | `?findOverlaps` |
| IR-COUNT-1 | `countOverlaps` = number of hits per query, same args as findOverlaps | `?countOverlaps` |
| IR-COUNT-2 | `countOverlaps` with empty subject is all zero; no subject = self | `?countOverlaps` |
| IR-COUNT-3 | R names the counts by query names; Python does not (ADR-0009) | `?countOverlaps` |
| IR-OVANY-1 | `overlapsAny` = count > 0 | `?overlapsAny` |
| IR-OVANY-2 | `overlapsAny` with zero-length query is zero-length | `?overlapsAny` |
| IR-SUBSET-1 | `subsetByOverlaps` keeps queries with a hit, preserving order, names, mcols | `?subsetByOverlaps` |
| IR-SUBSET-2 | `subsetByOverlaps(invert=TRUE)` keeps queries without a hit | `?subsetByOverlaps` |
