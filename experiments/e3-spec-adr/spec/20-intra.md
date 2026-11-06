# 20 Intra-range operations (INTRA)

Area: per-range transforms. All keep names and mcols (CONV.8), keep length.
Per-range arguments follow CONV.6. R-verified.

**INTRA.1** `shift(shift)`: `s' = s + shift`, width unchanged. Zero-width fine.
Ex: `shift({[1,3],[5,4]}, 2)` → `{[3,5],[7,6]}`.

**INTRA.2** `narrow(start=None, end=None, width=None)` selects a sub-range using positions *relative to the range*: relative position 1 is the range's first position, -1 its last (`-k` means `w-k+1`). Allowed argument sets: `start` alone (end kept), `end` alone (start kept), or any two of the three. `width` alone, or all three, is refused.
- Let `a`, `b` be the translated relative `start`, `end`. `start`+`end`: `s' = s+a-1`, `e' = s+b-1`. `start`+`width`: `s' = s+a-1`, `e' = s'+width-1`. `end`+`width`: `e' = s+b-1`, `s' = e'-width+1`. `start` alone: `e' = e`; `end` alone: `s' = s`.
- The result must satisfy `s' ≥ s`, `e' ≤ e`, `e' ≥ s'-1`, else refused. So a `start` translating to 0 is refused, `end=0` gives the zero-width `[s,s-1]`, and `start=w+1` gives the zero-width `[e+1,e]`.
Ex: `narrow({[10,13]}, start=2)` → `{[11,13]}`; `start=-2` → `{[12,13]}`; `end=-2` → `{[10,12]}`; `end=0` → `{[10,9]}`; `start=5` → `{[14,13]}`; `start=-5` refused; `start=4,end=2` refused; `narrow({[1,3]}, start=3, width=0)` → `{[3,2]}`.

**INTRA.3** `resize(width, fix="start")`; `width ≥ 0` else refused.
- `fix="start"`: `s' = s`.
- `fix="end"`: `e' = e` (`s' = e - width + 1`).
- `fix="center"`: `s' = s + floor((w - width)/2)` (floor toward -∞).
Ex: `resize({[1,5]}, 2, "center")` → `{[2,3]}`; `{[1,4]}`→3 → `{[1,3]}`; `{[1,3]}`→5 → `{[0,4]}`; `{[1,3]}`→0 end → `{[4,3]}`.

**INTRA.4** `flank(width, start=True, both=False)`. `start` may be per-range bool. Let `k = width`.
- `both=False`, `start=True`: `k ≥ 0` → `[s-k, s-1]`; `k < 0` → `[s, s-k-1]`.
- `both=False`, `start=False`: `k ≥ 0` → `[e+1, e+k]`; `k < 0` → `[e+k+1, e]`.
- `both=True`: uses `|k|`; `start=True` → `[s-|k|, s+|k|-1]`; `start=False` → `[e-|k|+1, e+|k|]`.
Zero-width input uses its own `s`, `e = s-1`.
Ex: `flank({[1,3]}, 2)` → `{[-1,0]}`; `flank({[1,3]}, 2, start=False)` → `{[4,5]}`; `flank({[1,3]}, -2)` → `{[1,2]}`; `flank({[1,3]}, 2, both=True)` → `{[-1,2]}`.

**INTRA.5** `promoters(upstream=2000, downstream=200)`: `[s-upstream, s+downstream-1]`; negative arguments refused.
Ex: `promoters({[1,3]}, 3, 1)` → `{[-2,1]}`.
