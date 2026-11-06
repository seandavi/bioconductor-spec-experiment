# 60 Overlaps and Hits (OVL)

Area: `find_overlaps` and friends; the `Hits` result. `q` = query range,
`s` = subject range. R-verified.

**OVL.1** Definitions. `gap(q,s)` per CONV notation. `ow(q,s) = max(0, min(e_q,e_s) - max(s_q,s_s) + 1)` (overlap width).

**OVL.2** `find_overlaps(subject, max_gap=-1, min_overlap=0, type="any", select="all")`. `max_gap ≥ -1` and `min_overlap ≥ 0`, else refused. A pair `(q,s)` is a hit iff the type condition (OVL.3–OVL.7) holds and `ow(q,s) ≥ min_overlap`.

**OVL.3** `type="any"`: hit iff `gap(q,s) ≤ max_gap`. With the default `max_gap=-1` this is `s_q ≤ e_s and s_s ≤ e_q`. Refused when `max_gap ≥ 0` and `min_overlap > 0`.
Consequences: a zero-width range overlaps a range strictly containing its position (`[5,4]` hits `[1,10]`, hits `[4,5]`) but not a range that starts or ends exactly at the position (`[5,4]` vs `[5,6]` and `[3,4]`: no hit); two zero-width ranges never hit at `max_gap=-1`. `max_gap=0` also accepts adjacent pairs; `min_overlap ≥ 1` excludes any pair involving a zero-width range.

**OVL.4** For types other than `any`, let `g = max(max_gap, 0)`.
- `start`: `|s_q - s_s| ≤ g`. `end`: `|e_q - e_s| ≤ g`. `equal`: both.
- No overlap is required (adjacent or disjoint pairs can hit; `min_overlap` can add it).
Ex: `start`, `max_gap=1`: `[4,5]` hits `[3,3]`, `[5,7]`; `[0,0]` hits `[1,10]`.

**OVL.5** `type="within"`: hit iff `s_s ≤ s_q and e_q ≤ e_s` (q inside s), and, when `max_gap ≥ 1`, also `w_s - w_q ≤ max_gap`. For `max_gap ≤ 0` the width difference is unrestricted (R quirk). ADR-0008.
Ex: q=`[5,8]`, subjects `[4,9],[3,10],[5,8],[3,9]`: `max_gap=0` → all four; `1` → `[5,8]`; `2` → `[4,9]`,`[5,8]`; `3` → `[4,9]`,`[5,8]`,`[3,9]`.

**OVL.6** `min_overlap` applies to every type as an extra filter (`ow ≥ min_overlap`; `0` = no filter).
Ex: q=`[5,8]`, subject `[5,9]`, `type="start"`: hit at `min_overlap=4`, none at `5`.

**OVL.7** `Hits` order: pairs ascending by `(query index, subject index)`. ADR-0009 (R orders subjects within a query by interval-index traversal; not replicated).

**OVL.8** `select="first"`/`"last"`: per query the smallest / largest subject index among its hits, NA if none. `select="arbitrary"`: any one hit's subject index (this implementation: the smallest), NA if none. Results are masked int64 of length `len(query)` (CONV.4). `select="all"` returns `Hits`.

**OVL.9** `Hits` holds `query_hits`, `subject_hits` (int64, equal length), `query_length`, `subject_length`; `len(hits)` = number of pairs.

**OVL.10** `count_overlaps(subject, ...)`: per query, the number of hits (same arguments as OVL.2). `overlaps_any(subject, ...)`: per query bool, hit count `> 0`. `subset_by_overlaps(subject, ..., invert=False)`: the queries with at least one hit, original order and duplicates kept; `invert=True` the others.

**OVL.11** Empty query or subject → no hits (counts 0, selects NA).
