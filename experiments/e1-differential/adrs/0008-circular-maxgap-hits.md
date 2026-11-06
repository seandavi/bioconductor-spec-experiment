# ADR-0008: circular sequences with maxgap > 0 on start / end / equal / within

Status: accepted (divergence from R)

## Context
On circular sequences R treats coordinates as periodic. For `type="any"` (all `maxgap`/`minoverlap`) and for
`start`/`end`/`equal`/`within` with the default `maxgap`, Python reproduces R exactly (hundreds of generated
cases): a query hits a subject if it hits any copy `subject + k*L`, ties and order as probed.
With `maxgap > 0` on the other types, R's answer depends on the internal engine path: the same
(query, subject) pair hits when searched alone, or among a different number of ranges, and misses in
another context (`div-circ-end-maxgap-*`: one pair alone misses; against five subjects the hit set is
{1,3,4,5}; with five identical queries it is {2,3,5}). No single rule over the pair reproduces this.

## Options
1. Reverse-engineer the two engines' circular wrap handling for every context.
2. Define one consistent rule: hit iff any copy satisfies the linear predicate.

## Decision
Option 2. A pairwise predicate that does not depend on how many other ranges are in the call is the
only behaviour users can reason about. The three `div-circ-end-maxgap-*` cases pin Python's answers;
the generators avoid this combination elsewhere.
