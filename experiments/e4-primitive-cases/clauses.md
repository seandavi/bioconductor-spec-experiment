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
