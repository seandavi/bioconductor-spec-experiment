# 40 Set operations (SET)

Area: `union`, `intersect`, `setdiff` and parallel `punion`, `pintersect`,
`psetdiff`. R-verified.

**SET.1** `union(y)`, `intersect(y)`, `setdiff(y)` treat each input as the set of positions it covers (zero-width ranges cover nothing). They return the set `x∪y`, `x∩y`, `x∖y` as maximal runs of positions (adjacent runs merge), ascending, never zero-width. No names/mcols.
Ex: `a={[1,5],[10,15],[20,25]}`, `b={[4,11],[12,13],[30,35]}`: union → `{[1,15],[20,25],[30,35]}`; intersect → `{[4,5],[10,13]}`; setdiff(a,b) → `{[1,3],[14,15],[20,25]}`.
Ex: `union({[1,3]},{[4,6]})` → `{[1,6]}`; `intersect({[1,3]},{[4,6]})` → `{}`; `union({[5,4]},{[5,4]})` → `{}`; `intersect({[1,10]},{[5,4],[8,12]})` → `{[8,10]}`.

**SET.2** Parallel versions require `len(x) == len(y)` (else refused) and work pairwise; names come from `x`; result length is `n`.

**SET.3** `punion(y, fill_gap=False)`: pair result is the hull `[min s, max e]`. If `gap(x_i,y_i) > 0` the pair is refused unless `fill_gap=True`.
Ex: `{[1,3]}∪{[4,6]}` → `{[1,6]}`; `{[1,3]}∪{[7,9]}` refused, with `fill_gap` → `{[1,9]}`; `{[5,4]}∪{[2,4]}` → `{[2,4]}`; `{[5,4]}∪{[1,3]}` refused (gap 1).

**SET.4** `pintersect(y, resolve_empty="none")`: let `s = max(s_x,s_y)`, `e = min(e_x,e_y)`. If `e ≥ s-1` the pair result is `[s,e]` (zero-width when the pair is adjacent or one is zero-width). If `e < s-1` (a real gap): `"none"` refuses; `"max.start"` gives `[s,s-1]`; `"start.x"` gives `[s_x,s_x-1]`.
Ex: `{[1,3]}∩{[4,6]}` → `{[4,3]}`; `{[1,10]}∩{[5,4]}` → `{[5,4]}`; `{[1,3]}∩{[7,9]}` refused; `max.start` → `{[7,6]}`; `start.x` → `{[1,0]}`.

**SET.5** `psetdiff(y)`: pair result is the positions of `x_i` not in `y_i`.
- One run `[a,b]` → `[a,b]`. No positions left (y covers x, or x is zero-width) → zero-width at `s_x`: `[s_x,s_x-1]`. y disjoint from x → `x_i` unchanged.
- Two runs (y strictly inside x) → refused. A zero-width y with `s_x < s_y ≤ e_x` also splits x and is refused; zero-width y with `s_y ≤ s_x` or `s_y > e_x` leaves x unchanged.
Ex: `{[1,10]}∖{[1,4]}` → `{[5,10]}`; `{[1,10]}∖{[8,20]}` → `{[1,7]}`; `{[4,5]}∖{[1,10]}` → `{[4,3]}`; `{[1,10]}∖{[4,5]}` refused; `{[1,10]}∖{[5,4]}` refused; `{[1,10]}∖{[1,0]}` → `{[1,10]}`.
