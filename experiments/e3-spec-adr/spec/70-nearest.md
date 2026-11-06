# 70 Nearest, precede, follow, distance (NEAR)

Area: neighbor queries between ranges. R-verified.

**NEAR.1** `distance(y)` is pairwise: `max(0, gap(x_i,y_i))`. Lengths equal, or one side length 1 (recycled); otherwise refused. Overlapping and adjacent ranges both give 0.
Ex: `[1,3]`–`[4,6]` → 0; `[1,3]`–`[5,6]` → 1; `[1,3]`–`[2,5]` → 0; `[5,4]`–`[1,3]` → 1; `[11,10]`–`[1,10]` → 0; `[12,11]`–`[1,10]` → 1.

**NEAR.2** `nearest(subject=None, select="arbitrary")`: for each query, the candidate set is every subject at minimal `distance` (adjacent and overlapping tie at 0; zero-width subjects count). `select="all"` returns `Hits` of every candidate (OVL.7 order); `"arbitrary"` returns one candidate index (this implementation: lowest). No subjects → NA. `subject=None` means the subject is `x` itself with the query's own index excluded (duplicates of it are not excluded; a length-1 `x` gives NA). Results for zero-width `x` ranges with `subject=None` are implementation-defined (R returns the range itself, twice).
Ex: `S={[1,3],[5,7],[10,12],[20,25],[20,22]}`: query `[4,4]` → candidates 0,1 (both adjacent); `[8,9]` → 1,2; `[15,15]` → 2 (gap 2 vs 4); `[21,21]` → 3,4 (overlap both); `[5,4]` → 1.

**NEAR.3** `precede(subject=None, select="first")`: candidates are subjects with `s_s > e_q`; the nearest-by-start group is the one with the smallest `s_s`. `"first"` returns the lowest index in that group, `"all"` returns `Hits` of the whole group. None → NA. Zero-width queries and subjects take part with their own `s`/`e`. `subject=None` → `x` itself.
Ex: `S` as NEAR.2: query `[4,4]` → 1; `[13,14]` → 3 (group 3,4); `[26,30]` → NA; `[5,4]` → 1; `[1,3]` vs subject `[3,6]` → NA.

**NEAR.4** `follow(subject=None, select="first")`: candidates are subjects with `e_s < s_q`; the group is the one with the largest `e_s`. `"first"` returns the *highest* index in the group (R quirk), `"all"` the whole group. ADR-0010.
Ex: query `[10,12]`, `S={[1,3],[3,5],[7,9],[8,9]}` → group 2,3 → first=3. Query `[26,30]` vs `S` of NEAR.2 → 3.
