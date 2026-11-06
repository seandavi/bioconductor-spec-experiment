# Red team: edges

46 cases: 46 agree, 0 strict expected failures. Accepted divergences are pinned to an ADR.

Legend: coordinates are start/width; hits/indices are R's 1-based values; NA is missing.
Tables abbreviate outputs; the comparator checks complete metadata, strand, and seqinfo.
'Agrees' includes matching errors. 'ADR judgment' is unresolved, not permission to diverge.

| Case id | R | Python | View |
|---|---|---|---|
| `rt-edge-ir-default-empty` | empty | same | Agrees |
| `rt-edge-gr-default-empty` | empty; levels=[] | same | Agrees |
| `rt-edge-ir-empty-reduce` | empty; revmap=[] | same | Agrees |
| `rt-edge-ir-empty-disjoin` | empty; revmap=[] | same | Agrees |
| `rt-edge-ir-empty-range` | empty | same | Agrees |
| `rt-edge-ir-empty-coverage` | runs=[]; lengths=[] | same | Agrees |
| `rt-edge-ir-empty-sort` | empty; names=[]; score=[] | same | Agrees |
| `rt-edge-ir-empty-unique` | empty; names=[]; score=[] | same | Agrees |
| `rt-edge-ir-empty-nearest` | no hits | same | Agrees |
| `rt-edge-ir-empty-union` | empty | same | Agrees |
| `rt-edge-gr-empty-reduce` | empty; revmap=[]; levels=[chrZ,chrA] | same | Agrees |
| `rt-edge-gr-empty-disjoin` | empty; revmap=[]; levels=[chrZ,chrA] | same | Agrees |
| `rt-edge-gr-empty-range` | empty; levels=[chrZ,chrA] | same | Agrees |
| `rt-edge-gr-empty-coverage` | chrZ: runs=[]; lengths=[]; chrA: runs=[0]; lengths=[4] | same | Agrees |
| `rt-edge-gr-empty-sort` | empty; names=[]; score=[]; levels=[chrZ,chrA] | same | Agrees |
| `rt-edge-gr-empty-unique` | empty; names=[]; score=[]; levels=[chrZ,chrA] | same | Agrees |
| `rt-edge-gr-empty-nearest` | no hits | same | Agrees |
| `rt-edge-gr-empty-union` | empty; levels=[chrZ,chrA] | same | Agrees |
| `rt-edge-zero-seqlength-trim` | 1/0,1/0,1/0,1/0,4/0; names=[a,b,c,d,e]; score=[1,2,3,4,5]; strand=[+,-,*,+,*]; levels=[chrZ,chrA,unused] | same | Agrees |
| `rt-edge-zero-seqlength-gaps` | 1/3,1/3,1/3,1/2,1/2,1/2; strand=[+,-,*,+,-,*]; levels=[chrZ,chrA,unused] | same | Agrees |
| `rt-edge-zero-seqlength-coverage` | chrZ: runs=[]; lengths=[]; chrA: runs=[0]; lengths=[3]; unused: runs=[0]; lengths=[2] | same | Agrees |
| `rt-edge-ordering-order-ascending` | [4,5,3,2,1,6] | same | Agrees |
| `rt-edge-ordering-order-descending` | [1,6,2,3,4,5] | same | Agrees |
| `rt-edge-ordering-sort-descending` | 5/0,5/0,5/0,5/0,5/0,5/0; names=[a,f,b,c,d,e]; score=[6,1,5,4,3,2]; strand=[*,*,*,-,+,+]; levels=[chrZ,chrA,unused] | same | Agrees |
| `rt-edge-ordering-sort-ignore-strand` | 5/0,5/0,5/0,5/0,5/0,5/0; names=[b,c,d,e,a,f]; score=[5,4,3,2,6,1]; strand=[*,-,+,+,*,*]; levels=[chrZ,chrA,unused] | same | Agrees |
| `rt-edge-ordering-rank-last` | [6,4,3,2,1,5] | same | Agrees |
| `rt-edge-ordering-duplicated-from-last` | [true,false,false,true,false,false] | same | Agrees |
| `rt-edge-ordering-unique-from-last` | 5/0,5/0,5/0,5/0; names=[b,c,e,f]; score=[5,4,2,1]; strand=[*,-,+,*]; levels=[chrZ,chrA,unused] | same | Agrees |
| `rt-edge-ordering-subset-empty-name` | error: subscript contains invalid names | also errors | Agrees |
| `rt-edge-ir-negative-nearest` | 1→1, 2→3, 3→3, 4→4 | same | Agrees |
| `rt-edge-ir-negative-precede` | 1→2, 2→4, 3→4, 4→4 | same | Agrees |
| `rt-edge-ir-negative-follow` | 1→1, 2→2, 3→3, 4→3 | same | Agrees |
| `rt-edge-ir-coverage-typed-na` | runs=[1,NA]; lengths=[1,5] | same | Agrees |
| `rt-edge-ir-coverage-overlap-na` | runs=[NA]; lengths=[6] | same | Agrees |
| `rt-edge-gr-negative-nearest` | 1→3, 2→3, 4→3 | 1→4, 2→3, 3→2, 4→3 | Accepted divergence (ADR-0007) |
| `rt-edge-gr-negative-precede` | no hits | 1→4, 3→2 | Accepted divergence (ADR-0007) |
| `rt-edge-gr-negative-follow` | error: logical subscript contains NAs | 2→1, 4→3 | Accepted divergence (ADR-0007) |
| `rt-edge-gr-coverage-typed-na` | chr1: runs=[1,NA]; lengths=[1,5] | same | Agrees |
| `rt-edge-gr-coverage-overlap-na` | chr1: runs=[NA]; lengths=[6] | same | Agrees |
| `rt-edge-trim-circularity-unknown` | 1/3; strand=[*] | same | Agrees |
| `rt-edge-trim-circularity-linear` | 1/3; strand=[*] | same | Agrees |
| `rt-edge-trim-circularity-circular` | -2/8; strand=[*] | same | Agrees |
| `rt-edge-all-zero-range` | -2/10 | same | Agrees |
| `rt-edge-all-zero-gaps` | -5/16 | same | Agrees |
| `rt-edge-all-zero-reduce` | -2/0,4/0,8/0; revmap=[[2], [4], [1, 3]] | same | Agrees |
| `rt-edge-all-zero-disjoin` | empty; revmap=[] | same | Agrees |
