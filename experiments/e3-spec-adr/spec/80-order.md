# 80 Ordering (ORD)

Area: `order`, `sort`, `rank`, `duplicated`, `unique` for `IRanges`
(GRanges: GR.15). R-verified.

**ORD.1** Ranges compare by `(start, end)` ascending (equivalently `(start, width)`). Two ranges are *equal* iff same start and width; names and mcols are ignored.

**ORD.2** `order(decreasing=False)`: indices (0-based) that sort ascending by ORD.1; ties keep index order. `decreasing=True` sorts descending by ORD.1, ties still in ascending index order.
Ex: `x={[5,9],[1,4],[5,8],[3,3],[1,4],[5,9],[5,4]}`: order → `(1,4,3,6,2,0,5)`; decreasing → `(0,5,2,6,3,1,4)`.

**ORD.3** `sort(decreasing=False)` = `x[order(decreasing)]`; names and mcols follow.

**ORD.4** `rank(ties_method="average")`: the rank value (1-based) of each range's key among all keys. `average`: mean of the tied ranks (float); `min`, `max`, `first` (ties by index): ints.
Ex: same `x` → average `(6.5,1.5,5,3,1.5,6.5,4)`; first `(6,1,5,3,2,7,4)`; min `(6,1,5,3,1,6,4)`.

**ORD.5** `duplicated(from_last=False)`: bool per range, true iff an equal range occurs earlier (later when `from_last`). `unique()` keeps the non-duplicated ranges in their original order (not sorted) with their names and mcols.
Ex: `x` above → duplicated `(F,F,F,F,T,T,F)`; `from_last` `(T,T,F,F,F,F,F)`; unique → `{[5,9],[1,4],[5,8],[3,3],[5,4]}`.
