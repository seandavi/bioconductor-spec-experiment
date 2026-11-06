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
| IR-NEAREST-1 | `nearest` (default `select="arbitrary"`) returns one subject among the minimum-distance ones (ADR-0005) | `?nearest` |
| IR-NEAREST-2 | `nearest(select="all")` returns every subject at the minimum distance, sorted by query then subject | `?nearest` |
| IR-NEAREST-3 | distance is 0 for overlapping and for adjacent subjects, so they tie | `?nearest` |
| IR-NEAREST-4 | equidistant subjects on both sides all tie | `?nearest` |
| IR-NEAREST-5 | zero-width queries/subjects use `[s, s-1]` coordinates | `?nearest` |
| IR-NEAREST-6 | empty subject gives NA (select arbitrary) or no hits (all); empty query gives zero-length | `?nearest` |
| IR-NEAREST-7 | `nearest` accepts only `select` arbitrary/all; others are errors | `?nearest` |
| IR-PRECEDE-1 | `precede` (default `select="first"`): subject with smallest start > end(x); lowest index among ties; NA if none | `?precede` |
| IR-PRECEDE-2 | `precede(select="all")` returns all subjects tied on that start | `?precede` |
| IR-PRECEDE-3 | overlapping subjects are excluded; adjacent ones count (start = end(x)+1) | `?precede` |
| IR-PRECEDE-4 | empty subject gives NA; `select` other than first/all is an error | `?precede` |
| IR-FOLLOW-1 | `follow` (default `select="last"`): subject with largest end < start(x); highest index among ties; NA if none | `?follow` |
| IR-FOLLOW-2 | `follow(select="all")` returns all subjects tied on that end | `?follow` |
| IR-FOLLOW-3 | ties on end; adjacent subjects count (end = start(x)-1) | `?follow` |
| IR-FOLLOW-4 | empty subject gives NA; `select` other than last/all is an error | `?follow` |
| IR-DISTANCE-1 | `distance` is `max(0, max(starts) - min(ends) - 1)` per pair (0 if overlapping or adjacent) | `?distance` |
| IR-DISTANCE-2 | length-1 argument recycles | `?distance` |
| IR-DISTANCE-3 | zero-length inputs give zero-length output | `?distance` |
| IR-DISTANCE-4 | R recycles cyclically with a warning; Python errors (ADR-0007) | `?distance` |
| IR-COV-1 | `coverage` is per-position depth for positions 1..max(end), run-length encoded (adjacent equal runs merge) | `?coverage` |
| IR-COV-2 | zero-width ranges add no depth, but their end still counts toward the result length | `?coverage` |
| IR-COV-3 | coverage of zero-length input is zero-length (or `width` zeros) | `?coverage` |
| IR-COV-4 | `shift` moves ranges before counting (scalar or per-range); length grows accordingly | `?coverage` |
| IR-COV-5 | positions below 1 are dropped | `?coverage` |
| IR-COV-6 | `width` pads with zeros or truncates; `width=0` gives zero-length | `?coverage` |
| IR-COV-7 | `weight` multiplies each range's contribution (scalar or per-range; negative and fractional allowed) | `?coverage` |
| IR-COV-8 | R recycles weight cyclically with a warning; Python errors (ADR-0007) | `?coverage` |
| IR-ORDER-1 | `order` sorts by start, then width | `?"IRanges-class"` |
| IR-ORDER-2 | `order(decreasing=TRUE)` reverses both keys | `?"IRanges-class"` |
| IR-ORDER-3 | ties keep input order, also when decreasing | `?"IRanges-class"` |
| IR-ORDER-4 | `order` of zero-length is zero-length | `?"IRanges-class"` |
| IR-SORT-1 | `sort` reorders by `order`, carrying names and mcols | `?"IRanges-class"` |
| IR-SORT-2 | `sort` of zero-length is zero-length | `?"IRanges-class"` |
| IR-RANK-1 | `rank` by (start, width) with `ties.method` average/first/last/min/max; `last` gives later ties lower ranks | `?"IRanges-class"` |
| IR-RANK-2 | `rank` of zero-length is zero-length | `?"IRanges-class"` |
| IR-RANK-3 | `ties.method` has no default: omitting it, or an unknown value, is an error | `?"IRanges-class"` |
| IR-DUP-1 | `duplicated` marks later repeats of the same (start, width) | `?"IRanges-class"` |
| IR-DUP-2 | `duplicated(fromLast=TRUE)` marks earlier repeats | `?"IRanges-class"` |
| IR-DUP-3 | `duplicated` ignores names and mcols; zero-width ranges compare by start | `?"IRanges-class"` |
| IR-DUP-4 | `duplicated` of zero-length is zero-length | `?"IRanges-class"` |
| IR-UNIQUE-1 | `unique` keeps first occurrences with their names and mcols | `?"IRanges-class"` |
| IR-UNIQUE-2 | `unique` of zero-length is zero-length | `?"IRanges-class"` |
| GR-CON-1 | `GRanges(seqnames, ranges, strand)` stores seqnames, ranges and strand per range | `?GRanges` |
| GR-CON-2 | strand defaults to `*`; length-1 seqnames/strand recycle | `?GRanges` |
| GR-CON-3 | without `seqinfo`, seqlevels are the seqnames in order of first appearance | `?GRanges` |
| GR-CON-4 | names and mcols are carried through; ranges may be given by start/end | `?GRanges` |
| GR-CON-5 | zero-length and zero-width GRanges are valid | `?GRanges` |
| GR-CON-6 | with `seqinfo`, seqlevels/lengths/circularity/genome come from it (unused levels kept) | `?GRanges` |
| GR-CON-7 | seqname absent from seqinfo, NA seqname, bad strand, length mismatch, negative width are errors | `?GRanges` |
| GR-ACC-1 | accessors: seqnames, strand, start, end, width, names, ranges, seqlevels, seqlengths, isCircular, genome, seqinfo, length | `?GRanges` |
| SI-CON-1 | `Seqinfo(seqnames, seqlengths, isCircular, genome)`; omitted fields are NA | `?Seqinfo` |
| SI-CON-2 | a scalar `genome` recycles to every sequence | `?Seqinfo` |
| SI-CON-3 | duplicate names, negative length, vector length mismatch are errors | `?Seqinfo` |
| SI-ACC-1 | `seqlevels`, `seqlengths`, `isCircular`, `genome` on Seqinfo (named vectors) | `?Seqinfo` |
| GR-SUB-1 | integer/logical subscripts select ranges; repeats allowed | `?GRanges` |
| GR-SUB-2 | character subscripts select by names | `?GRanges` |
| GR-SUB-3 | subsetting keeps the full seqinfo (including unused seqlevels) | `?GRanges` |
| GR-SUB-4 | out-of-bounds index or unknown name is an error | `?GRanges` |
| GR-CAT-1 | `c()` merges seqinfo: x's seqlevels first, then new ones from later arguments | `?"c,GenomicRanges-method"` |
| GR-CAT-2 | NA length/circularity/genome are filled from the other side | `?"c,GenomicRanges-method"` |
| GR-CAT-3 | conflicting non-NA seqlength, circularity or genome for the same sequence is an error | `?"c,GenomicRanges-method"` |
| GR-CAT-4 | no common seqlevels: warning only; levels are unioned | `?"c,GenomicRanges-method"` |
| GR-CAT-5 | `c()` pads names with `""` and fills missing mcols with NA | `?"c,GenomicRanges-method"` |
| GR-TRIM-1 | `trim` clamps ranges to `[1, seqlength]` on sequences of known length | `?"trim,GenomicRanges-method"` |
| GR-TRIM-2 | circular sequences and sequences of unknown length are not trimmed | `?"trim,GenomicRanges-method"` |
| GR-TRIM-3 | ranges wholly outside collapse to a zero-width range at the boundary (`[1,0]` or `[len+1,len]`) | `?"trim,GenomicRanges-method"` |
| GR-TRIM-4 | `trim` of zero-length is zero-length | `?"trim,GenomicRanges-method"` |
| GR-SHIFT-1 | `shift` moves ranges regardless of strand; keeps seqnames, strand, names, mcols, seqinfo | `?"intra-range-methods"` |
| GR-SHIFT-2 | `narrow` is strand-agnostic (coordinates relative to range start) | `?narrow` |
| GR-SHIFT-3 | shifting past a sequence end is allowed (no error, no trim) | `?"intra-range-methods"` |
| GR-RESIZE-1 | `resize` fix `start`/`end` refer to the 5'/3' end: swapped on `-` strand; `*` acts as `+` | `?resize` |
| GR-RESIZE-2 | `resize(fix="center")` is strand-independent | `?resize` |
| GR-RESIZE-3 | `ignore.strand=TRUE` treats all ranges as `+` | `?resize` |
| GR-RESIZE-4 | `use.names=FALSE` drops names; negative width is an error | `?resize` |
| GR-FLANK-1 | `flank(width>0)`: upstream of start for `+`/`*`; downstream of end for `-` (and swapped by `start=FALSE`) | `?flank` |
| GR-FLANK-2 | `flank(both=TRUE)` spans `abs(width)` each side of the strand-aware boundary | `?flank` |
| GR-FLANK-3 | `flank(width<0)` is the strand-aware region inside the range | `?flank` |
| GR-FLANK-4 | `ignore.strand=TRUE` treats all ranges as `+` | `?flank` |
| GR-FLANK-5 | flank may extend before position 1 (out of bounds, warning only) | `?flank` |
| GR-PROMOTERS-1 | `promoters`: `+`/`*` -> `[start-up, start+down-1]`; `-` -> `[end-down+1, end+up]` | `?promoters` |
| GR-PROMOTERS-2 | `promoters` defaults 2000/200 (may go out of bounds, warning only) | `?promoters` |
| GR-PROMOTERS-3 | negative upstream/downstream is an error | `?promoters` |
| GR-INTRA-1 | intra-range ops on zero-length GRanges return zero-length GRanges with seqinfo | `?"intra-range-methods"` |
| GR-RANGE-1 | `range` is computed per (seqlevel, strand) group; output sorted by seqlevel, strand (+ - *), position | `?range` |
| GR-RANGE-2 | `range(ignore.strand=TRUE)` merges strands per seqlevel; result strand `*` | `?range` |
| GR-RANGE-3 | `range` of zero-length is zero-length with seqinfo | `?range` |
| GR-REDUCE-1 | `reduce` runs per (seqlevel, strand); `*` is its own group (not merged with `+`/`-`) | `?reduce` |
| GR-REDUCE-2 | `reduce(ignore.strand=TRUE)` merges across strands; result strand `*` | `?reduce` |
| GR-REDUCE-3 | `min.gapwidth` and `drop.empty.ranges` apply per group as for IRanges | `?reduce` |
| GR-REDUCE-4 | `with.revmap=TRUE`: `revmap` holds indices into the whole input | `?reduce` |
| GR-REDUCE-5 | `reduce` of zero-length keeps seqinfo | `?reduce` |
| GR-REDUCE-6 | output order follows seqlevels order (not alphabetical), then strand + - * | `?reduce` |
| GR-DISJOIN-1 | `disjoin` per (seqlevel, strand) group | `?disjoin` |
| GR-DISJOIN-2 | `disjoin(ignore.strand=TRUE)` pools strands; result `*` | `?disjoin` |
| GR-DISJOIN-3 | `disjoin(with.revmap=TRUE)` indexes the whole input | `?disjoin` |
| GR-DISJOIN-4 | `disjoin` of zero-length keeps seqinfo | `?disjoin` |
| GR-GAPS-1 | `gaps` returns, for every seqlevel and each strand, the gaps within `[1, seqlength]` (unused sequences/strands give whole-sequence gaps) | `?gaps` |
| GR-GAPS-2 | `gaps(ignore.strand=TRUE)`: one `*` group per seqlevel | `?gaps` |
| GR-GAPS-3 | `start`/`end` are scalars applied to every group, overriding 1 / seqlength | `?gaps` |
| GR-GAPS-4 | zero-length input: whole-sequence gaps where lengths are known; unused level without length gives nothing | `?gaps` |
| GR-GAPS-5 | unknown seqlength: the end defaults to the group's coverage extent (as IRanges); groups without ranges yield nothing | `?gaps` |
| GR-GAPS-6 | ranges beyond seqlength are clipped; circular sequences are not special; seqlength 0 gives nothing | `?gaps` |
| GR-UNION-1 | `union` per (seqlevel, strand) group; strands `+`,`-`,`*` never mix | `?"setops"` |
| GR-UNION-2 | `union(ignore.strand=TRUE)` pools strands; result `*` | `?"setops"` |
| GR-UNION-3 | `*` is not a wildcard in set operations | `?"setops"` |
| GR-UNION-4 | set ops with a zero-length argument | `?"setops"` |
| GR-UNION-5 | result seqinfo is the merge of both; seqlevel order x first | `?"setops"` |
| GR-INTERSECT-1 | `intersect` per group | `?"setops"` |
| GR-INTERSECT-2 | `intersect(ignore.strand=TRUE)` | `?"setops"` |
| GR-INTERSECT-3 | `intersect` of `+` with `*` is empty | `?"setops"` |
| GR-INTERSECT-4 | `intersect` with a zero-length argument | `?"setops"` |
| GR-INTERSECT-5 | result seqinfo is the merge of both | `?"setops"` |
| GR-SETDIFF-1 | `setdiff` per group | `?"setops"` |
| GR-SETDIFF-2 | `setdiff(ignore.strand=TRUE)` | `?"setops"` |
| GR-SETDIFF-3 | `setdiff` of `+` minus `*` removes nothing | `?"setops"` |
| GR-SETDIFF-4 | `setdiff` with a zero-length argument | `?"setops"` |
| GR-SETDIFF-5 | result seqinfo is the merge of both | `?"setops"` |
| GR-PUNION-1 | `punion` pairs: compatible seqname and strand required; result keeps x's names/seqnames/strand, drops mcols | `?"setops"` |
| GR-PUNION-2 | `punion` result strand is x's even when x is `*` | `?"setops"` |
| GR-PUNION-3 | `punion` gap rule as IRanges (`fill.gap`) | `?"setops"` |
| GR-PUNION-4 | incompatible pair or length mismatch is an error | `?"setops"` |
| GR-PINTERSECT-1 | `pintersect` adds logical mcol `hit` (compatible pair with non-ambiguous intersection); non-hits are `[start(x), start(x)-1]`; x's mcols kept | `?"setops"` |
| GR-PINTERSECT-2 | adjacent and zero-width-inside pairs are hits with an empty range; zero-width outside is a miss | `?"setops"` |
| GR-PINTERSECT-3 | hit strand is x's, or y's when x is `*` | `?"setops"` |
| GR-PINTERSECT-4 | strand/seqname mismatch is a miss (`ignore.strand=TRUE` ignores strand only); no error | `?"setops"` |
| GR-PSETDIFF-1 | `psetdiff` subtracts compatible pairs; keeps x's names, drops mcols | `?"setops"` |
| GR-PSETDIFF-2 | incompatible pairs return x unchanged | `?"setops"` |
| GR-PSETDIFF-3 | compatible pair strand is x's, or y's when x is `*` | `?"setops"` |
| GR-PSETDIFF-4 | y strictly inside x is an error | `?"setops"` |
