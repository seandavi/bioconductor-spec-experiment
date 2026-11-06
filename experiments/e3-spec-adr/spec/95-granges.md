# 95 GRanges (GR)

Area: `GRanges` = `IRanges` + seqnames + strand + seqinfo. Everything
in 10–80 applies per range; this file states what changes. Strand codes
`+ < - < *` (order used for sorting). R-verified (GenomicRanges 1.65.1).

## Core

**GR.1** `GRanges(seqnames, ranges, strand=None, seqinfo=None, names=None, mcols=None)`. `seqnames`: one string per range. `ranges`: an `IRanges` of length `n` (a length-1 input is recycled to `n`; its names/mcols are ignored). `strand`: one of `+ - *` (scalar recycled, or per range); default `*`. `names` default: none. Lengths must otherwise agree, else refused.

**GR.2** `strand` values other than `+`, `-`, `*` are refused. Seqnames absent from an explicit `Seqinfo` are refused (SEQI.4).

**GR.3** Accessors: `seqnames` (strings, length `n`), `strand` (strings), `ranges` (bare `IRanges`, no names/mcols), `start`, `end`, `width`, `names`, `mcols`, `seqinfo`, `len`. Ranges outside `[1, seqlength]` are allowed at construction; no check is made (GR.17).

**GR.4** Subsetting follows CORE.9/CORE.10, carrying seqnames, strand, names, mcols; `seqinfo` is kept whole (SEQI.5).

**GR.5** `concat(a, b, ...)`: ranges, seqnames, strand appended; `seqinfo` merged left to right (SEQI.2); names and mcols per CORE.12. Seqnames are matched by name, not level index.
Ex: `concat(chr2:1-2:+, chr3:1-2:-)` → levels `(chr2,chr3)`.

**GR.6** Strand compatibility: two strands are *compatible* iff equal or either is `*`. Two ranges are *comparable* iff same seqname and compatible strands. `ignore_strand=True` makes strands compatible and (where an output strand is produced by grouping) outputs `*`.

## Intra-range (strand-aware)

**GR.7** `shift`, `narrow`: as INTRA.1/2, independent of strand.

**GR.8** `resize`, `flank`, `promoters` treat `*` like `+`; for `-` the roles of start and end swap (the 5' end of a `-` range is its `e`):
- `resize(fix="start")` on `-` keeps `e`; `fix="end"` keeps `s`; `"center"` is not strand-dependent.
- `flank(start=True)` on `-` flanks the `e` side (`[e+1,e+k]` for `k ≥ 0`); `start=False` flanks the `s` side; `both` and negative `k` follow the same mirror.
- `promoters` on `-`: `[e-downstream+1, e+upstream]`.
- `ignore_strand=True` treats every range as `+`.
Ex (`a:10-14`): `+`: `resize 2` → 10-11, `flank 3` → 7-9, `promoters(4,2)` → 6-11. `-`: `resize 2` → 13-14, `flank 3` → 15-17, `flank 3 start=False` → 7-9, `promoters(4,2)` → 13-18.

**GR.9** `trim()`: for a range on a level with known `seqlength L` that is not circular: `s' = min(max(s,1), L+1)`, `e' = max(min(e,L), 0)`. Other ranges are unchanged. Out-of-bound ranges can end zero-width.
Ex: L=10: `[12,15]` → `[11,10]`; `[-5,-3]` → `[1,0]`; `[-5,200]` → `[1,10]`; circular or unknown L: unchanged.

## Inter-range and set operations

**GR.10** Grouping. `range`, `reduce`, `disjoin` run the IRanges operation (INTER.*) separately on each group of ranges with the same `(seqname, strand)`; with `ignore_strand=True` the group is the seqname alone and the result strand is `*`. Output groups ascend by `(seqlevel position, strand code)` and ranges ascend within a group. `seqinfo` is kept; names/mcols are dropped. `revmap` indices refer to the input `GRanges` (0-based), listed per INTER.5/INTER.8. `range` emits one range per non-empty group.
Ex: `{b:1-3:+, a:5-6:+, a:7-9:+, a:5-6:-, a:6-8:*, b:2-6:+}` (levels b,a,c) reduce → `b:1-6:+, a:5-9:+, a:5-6:-, a:6-8:*`; `ignore_strand` → `b:1-6:*, a:5-9:*`.

**GR.11** `gaps(start=None, end=None, ignore_strand=False)`: for every level (level order) and every strand `+,-,*` (just `*` with `ignore_strand`) the gaps of that group's ranges within the window `[start, end]`. Defaults: `start = 1`; `end = seqlength` of the level. A group whose `end` is unknown (no `seqlength` and no `end` argument) contributes only if it has ranges, with window end = its largest `e` (so there is no trailing gap); empty groups contribute nothing. Groups with a known window and no ranges yield the whole window. Explicit `start`/`end` apply to all groups. Output order as GR.10.
Ex: levels `a:NA, b:30`, ranges `a:5-10:+, a:20-25:+, b:3-4:-` → `a:1-4:+, a:11-19:+, b:1-30:+, b:1-2:-, b:5-30:-, b:1-30:*`.

**GR.12** `union`, `intersect`, `setdiff` (SET.1) run per `(seqname, strand)` group; `*` is its own group, not a wildcard. Group result order as GR.10. `ignore_strand=True` pools strands and outputs `*`.
Ex: `A={a:1-8:+, a:10-20:-, b:1-5:+}`, `B={a:5-12:+, a:10-15:*, b:3-9:-}`: union → `a:1-12:+, a:10-20:-, a:10-15:*, b:1-5:+, b:3-9:-`; intersect → `a:5-8:+`; setdiff → `a:1-4:+, a:10-20:-, b:1-5:+`; union with `ignore_strand` → `a:1-20:*, b:1-9:*`.

**GR.13** Parallel `punion`, `pintersect`, `psetdiff` (SET.2–5) apply the IRanges rule to pairs whose members are comparable (GR.6); the result strand is `x`'s, unless `x` is `*`, then `y`'s (with `ignore_strand`, always `x`'s). For a non-comparable pair: `punion` refuses, `pintersect` yields the zero-width `[s_x,s_x-1]` on `x`'s seqname and strand, `psetdiff` yields `x_i` unchanged. With `ignore_strand`, only the seqname decides comparability.
Ex: `a:1-8:+` ∩ `a:5-12:*` → `a:5-8:+`; ∩ `a:5-12:-` → `a:1-0:+`; `a:1-8:*` ∩ `a:5-12:-` → `a:5-8:-`.

## Overlaps and neighbors

**GR.14** `find_overlaps` and its relatives (OVL.*) count a pair as a hit only if it is comparable (GR.6) *and* satisfies the IRanges rule. `ignore_strand` applies. `Hits` order as OVL.7. `max_gap` applies only within comparable pairs.
Ex: queries `a:5-10:{+,-,*}`, subjects `a:1-7:{+,-,*}` → `+`–`+`, `+`–`*`, `-`–`-`, `-`–`*`, `*`–all three.

**GR.15** `distance`: NA (masked) for non-comparable pairs, else NEAR.1. `nearest`, `precede`, `follow` choose among comparable subjects only. `precede`/`follow` for a `-` query mirror direction: `precede` uses subjects at lower coordinates (`e_s < s_q`, largest `e_s`), `follow` subjects at higher coordinates (`s_s > e_q`, smallest `s_s`); `*` queries (and `ignore_strand`) behave as `+`. Tie rules of NEAR.3/NEAR.4 stay attached to the operation (`precede first` = lowest index, `follow first` = highest). `nearest` is NEAR.2 over the comparable subjects (all ties kept; ADR-0011).
Ex: S=`a:1-3:+, a:10-12:+, a:20-22:-, a:30-32:*, b:5-7:+, a:40-42:-`; query `a:6-7:+` → precede 1, follow 0; `a:6-7:-` → precede NA, follow 2; `b:1-2:-` → nearest/precede/follow NA (only `b:+` exists: not comparable).

## Order and coverage

**GR.16** Key of a range: `(position of seqname in seqlevels, strand code, start, end)`, ascending; ties by index. `order`, `sort`, `rank`, `duplicated`, `unique` follow ORD.2–ORD.5 with this key; equality = same seqname, strand, start, width. `decreasing=True` reverses the key order (ties still ascending index). Seqlevels sort in level order, not alphabetically.
Ex: levels `(b,a,c)`; `{b:5-6:+, a:1-2:-, b:3-4:-, c:1-2:*, a:1-2:+, b:3-3:+, a:1-2:+}` → order `(5,0,2,4,6,1,3)`.

**GR.17** `coverage(shift=0, width=None)` ignores strand and returns a mapping seqlevel → `Rle`, one entry per level in level order. Length of a level: `width` if given; else its `seqlength` if known; else the largest `e+shift` among its ranges (0 if none). Ranges are truncated to that length (COV.1–COV.3 otherwise).
Ex: levels `a:15, b:5, c:4`, ranges `a:1-4, a:3-6, b:2-3, a:10-12` → `a: 2×1,2×2,2×1,3×0,3×1,3×0`; `b: 1×0,2×1,2×0`; `c: 4×0`.
