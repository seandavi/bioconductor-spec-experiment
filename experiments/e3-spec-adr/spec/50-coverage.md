# 50 Coverage (COV)

Area: `coverage` and its run-length result. R-verified.

**COV.1** `coverage(shift=0, width=None)`: for position `p ≥ 1`, the value is the number of ranges `r` with `s_r+shift ≤ p ≤ e_r+shift`. `shift` is a scalar or per-range (CONV.6). Positions `< 1` are discarded. Zero-width ranges cover nothing.

**COV.2** Length `L` of the result: `width` if given (ranges are truncated to `[1,width]`, short coverage is padded with 0); else `max(0, max(e_r+shift))` over all ranges, zero-width included; the empty input gives `L=0`.
Ex: `{[1,3],[5,4],[10,13],[10,10]}` → runs `3×1, 6×0, 1×2, 3×1` (L=13). `{[5,4]}` → `4×0`. `{[2,3],[4,6]}` with `width=10` → `1×0, 5×1, 4×0`; `width=4` → `1×0, 3×1`; `shift=1` → `2×0, 5×1`. `{[-2,3]}` → `3×1`.

**COV.3** The result is an `Rle`: parallel `lengths` and `values` (int64), runs maximal (adjacent runs differ in value), lengths `> 0`, `sum(lengths) = L`; `L=0` has no runs. `Rle` offers `len`, `to_numpy()` (decoded values) and equality on runs.
