# 90 Seqinfo (SEQI)

Area: sequence metadata shared by `GRanges`. R-verified.

**SEQI.1** A `Seqinfo` is an ordered list of distinct sequence names (seqlevels) with, per level: `seqlength` (int ≥ 0 or absent), `is_circular` (bool or absent), `genome` (str or absent). Level order is meaningful (GR.15).

**SEQI.2** Merge (`Seqinfo.merge(y)`, used by `concat`): levels = `x`'s levels in order, then `y`'s levels not in `x` in `y`'s order. For a level in both, each attribute: both present and different → refuse; one present → that value; both absent → absent. Levels present in only one side keep their attributes. No warning is required when the level sets differ.
Ex: `x=(a:10, b:20)`, `y=(c:30, a:absent)` → `(a:10, b:20, c:30)`; `x=(a:10)`, `y=(a:11)` refused; `x=(a, circular)`, `y=(a, circular absent)` → circular.

**SEQI.3** Accessors on a `GRanges`: `seqlevels`, `seqlengths`, `is_circular`, `genome` return sequences in level order, absent as `None`.
Ex: levels `(b,a,c)`, lengths `(100,200,None)`.

**SEQI.4** Default `Seqinfo` for a `GRanges` built from seqnames only: levels are the distinct seqnames in order of first appearance, all attributes absent. An explicit `Seqinfo` fixes the level order; a seqname not among its levels is refused.
Ex: seqnames `(chr2,chr1,chr10,chr1)` → levels `(chr2,chr1,chr10)`.

**SEQI.5** Levels with no ranges are kept by subsetting, sorting, and all operations; the `Seqinfo` of a result equals the input's unless stated (concat: SEQI.2).
