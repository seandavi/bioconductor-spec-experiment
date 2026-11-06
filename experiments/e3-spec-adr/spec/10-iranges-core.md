# 10 IRanges core (CORE)

Area: construction, accessors, names, mcols, subsetting, concatenation.
R-verified (IRanges 2.47.2).

## Construction

**CORE.1** An `IRanges` is built from exactly two of `start`, `end`, `width` (integer sequences). The third is derived: `end = start+width-1`; `width = end-start+1`; `start = end-width+1`.
Ex: `start=(1,5), end=(3,4)` → `{[1,3],[5,4]}`; `end=(3,4), width=(3,0)` → same.

**CORE.2** Refuse when a derived or given width is negative (`end < start-1`, `width < 0`).

**CORE.3** Inputs of length 1 are recycled to the longest input; other length mismatches are refused (R warns and recycles; ADR-0006).
Ex: `start=(1,2,3), width=2` → `{[1,2],[2,3],[3,4]}`.

**CORE.4** Missing (NA) start/end/width is refused.

**CORE.5** The empty `IRanges` (length 0) is valid.

## Accessors

**CORE.6** `start`, `end`, `width` return int64 arrays of length `n`; `len(x) = n`.

**CORE.7** `names` is absent (`None`) or one entry per range; an entry may be absent (`None`). Names need not be unique.

**CORE.8** `mcols` is a table with `n` rows and any number of named columns. No columns = no metadata. A table whose height is not `n` is refused. Replacing names or mcols returns a new object.

## Subsetting

**CORE.9** `x[i]` selects ranges by: an int (negative counts from the end); a slice; an int sequence (order and duplicates kept); a bool sequence of length `n`; a name; a sequence of names. The result is always an `IRanges` (an int yields length 1). Out-of-range or NA positions, wrong-length bool masks, and unknown names are refused.
Ex: `{[1,3],[5,4],[10,13]}[(2,0)]` → `{[10,13],[1,3]}`; `[-1]` → `{[10,13]}`.

**CORE.10** A name selects its first match. Subsetting carries names and mcols rows of the selected ranges. Subsetting by an int array that is empty yields the empty `IRanges`.

## Concatenation

**CORE.11** `concat(a, b, ...)` appends ranges in argument order. At least one object is required.

**CORE.12** Names: if no input has names the result has none; otherwise each element keeps its name and elements of unnamed inputs get `""`. mcols: columns are the union in order of first appearance; a row from an input lacking a column gets NA.
Ex: names `("a","b","c")` concat unnamed length-3 → `("a","b","c","","","")`.
