# 30 Inter-range operations (INTER)

Area: `range`, `reduce`, `gaps`, `disjoin`. Results carry no names/mcols
(CONV.8) except the `revmap` column. R-verified.

**INTER.1** `range()` returns one range `[min s, max e]` over all ranges (zero-width included); the empty input gives the empty `IRanges`.
Ex: `range({[5,4],[1,0]})` → `{[1,4]}`.

**INTER.2** `reduce(min_gapwidth=1)` merges ranges. Sort by `(s, e)` (ties by input index). Scan in that order keeping a current group `[gs, ge]`; a range joins the group iff `s ≤ ge + min_gapwidth`, then `ge = max(ge, e)`; otherwise it starts a new group. Output the groups in order. `min_gapwidth ≥ 0`, else refused.
- `min_gapwidth=1` merges overlapping and adjacent ranges; `0` merges only overlapping ranges (`s ≤ ge`).
Ex (default): `{[1,3],[4,6],[8,8],[9,12],[20,25],[5,5]}` → `{[1,6],[8,12],[20,25]}`. `min_gapwidth=0` on `{[1,3],[4,6]}` → unchanged. `min_gapwidth=3` on `{[1,6],[8,12],[20,25]}` → `{[1,12],[20,25]}`.

**INTER.3** Zero-width ranges follow INTER.2 unchanged: they can seed a group, join one, and survive alone.
Ex: `{[1,3],[5,4],[10,13],[10,10]}` → `{[1,3],[5,4],[10,13]}`. `{[1,10],[5,4],[5,5],[6,5]}` → `{[1,10]}`. `{[1,3],[4,3]}` → `{[1,3]}`; with `min_gapwidth=0` → `{[1,3],[4,3]}`.
Quirk: with `min_gapwidth=0`, a zero-width range equal to the preceding scanned range, and not absorbed into a group, vanishes (output and revmap): `{[14,13],[14,13],[14,20]}` → `{[14,13],[14,20]}` with revmap `((0),(2))`. ADR-0013.

**INTER.4** `reduce(drop_empty_ranges=True)`: a zero-width range that would *seed* a new group (no group in progress, or it does not join the current group) is discarded instead; a zero-width range that joins an existing group is processed as in INTER.3 (it can extend `ge`). Groups therefore never have zero width.
Ex: `{[1,3],[5,4],[10,13]}` → `{[1,3],[10,13]}`. `{[5,4],[6,8]}` with `min_gapwidth=2` → `{[6,8]}` (without the flag: `{[5,8]}`). `{[1,3],[5,4]}` with `min_gapwidth=3` → `{[1,4]}`.
ADR-0007.

**INTER.5** `reduce(with_revmap=True)` adds mcols column `revmap`: for each output range, the input indices (0-based) of the ranges in its group, in the INTER.2 scan order `(s, e, index)`. With `drop_empty_ranges=True`, zero-width ranges (discarded or joined) appear in no list, though a joined one can still extend its group's end.
Ex: `{[1,3],[4,6],[8,8],[9,12],[20,25],[5,5]}` → revmap `((0,1,5),(2,3),(4))`. `{[1,3],[5,4],[10,13],[10,10]}` → `((0),(1),(3,2))`; with `drop_empty_ranges`, `{[7,6],[14,13],[3,5],[19,19],[12,16]}` → `{[3,5],[12,16],[19,19]}`, revmap `((2),(4),(3))`.

**INTER.6** `gaps(start=None, end=None)` returns the maximal runs of positions in `[start,end]` covered by no range, ascending, never zero-width. Zero-width ranges cover nothing. Ranges outside the window are irrelevant. Defaults (when not given): `start = min s`, `end = max e` over the *positive-width* ranges of `x`. If `x` has none and a bound is missing, the result is empty.
Ex: `x={[3,5],[7,10],[15,18]}`: default → `{[6,6],[11,14]}`; `start=1,end=20` → `{[1,2],[6,6],[11,14],[19,20]}`; `start=20,end=25` → `{[20,25]}`; `start=8,end=9` → `{}`; `start=12` → `{[12,14]}`; `end=12` → `{[6,6],[11,12]}`. `gaps({}, start=1, end=5)` → `{[1,5]}`. `gaps({[3,5],[6,5],[8,10]}, 1, 12)` → `{[1,2],[6,7],[11,12]}`; `gaps({[1,2],[5,5],[13,12]}, start=1)` → `{[3,4]}` (the zero-width range does not set `end`).

**INTER.7** `disjoin(with_revmap=False)`: let `B` = sorted distinct values of every `s` and every `e+1` over *all* ranges (zero-width included). The output is every interval `[b_i, b_{i+1}-1]` for consecutive `b_i, b_{i+1}` in `B` that is covered by at least one positive-width range. Adjacent pieces are not merged. A zero-width range inside a covered stretch therefore splits it.
Ex: `{[1,5],[3,8],[10,12],[5,4],[20,19],[6,6]}` → `{[1,2],[3,4],[5,5],[6,6],[7,8],[10,12]}`. `{[1,3],[4,6]}` → unchanged. `{[1,10],[5,4]}` → `{[1,4],[5,10]}`. `{[5,4]}` → `{}`.

**INTER.8** `disjoin(with_revmap=True)` adds `revmap`: per output piece, the 0-based indices (ascending) of the positive-width input ranges covering it.
Ex: `{[1,5],[3,8],[10,12]}` → pieces `{[1,2],[3,5],[6,8],[10,12]}`, revmap `((0),(0,1),(1),(2))`.
