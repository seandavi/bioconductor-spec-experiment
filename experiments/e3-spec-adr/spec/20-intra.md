# 20 Intra-range operations (INTRA)

Area: per-range transforms. All keep names and mcols (CONV.8), keep length.
Per-range arguments follow CONV.6. R-verified.

**INTRA.1** `shift(shift)`: `s' = s + shift`, width unchanged. Zero-width fine.
Ex: `shift({[1,3],[5,4]}, 2)` → `{[3,5],[7,6]}`.

**INTRA.2** `narrow(start=None, end=None, width=None)` selects a sub-range using positions *relative to the range*: relative position 1 is the range's first position, -1 its last. Allowed argument sets: `start` alone (end kept), `end` alone (start kept), or any two of the three. `width` alone, or relative position 0, is refused.
- `start` and `end` give relative positions `a`, `b` (negatives translated: `-k` → `width-k+1`): `s' = s+a-1`, `e' = s+b-1`.
- `start`+`width`: `s' = s+a-1`, `w' = width`. `end`+`width`: `e' = s+b-1`, `w' = width`.
- The result must lie inside the original range, zero-width results included at the boundary: `a ≤ width+1`, `b ≤ width`, `w' ≥ 0`. Otherwise refused.
Ex: `narrow({[10,13]}, start=2)` → `{[11,13]}`; `start=-2` → `{[12,13]}`; `end=-2` → `{[10,12]}`; `narrow({[1,3]}, start=3, width=0)` → `{[3,2]}`.

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
