# Red team: circular

33 cases: 16 agree, 17 strict expected failures. No accepted differences.

Legend: coordinates are start/width; hits/indices are R's 1-based values; NA is missing.
Tables abbreviate outputs; the comparator checks complete metadata, strand, and seqinfo.
'Agrees' includes matching errors. 'ADR judgment' is unresolved, not permission to diverge.

| Case id | R | Python | View |
|---|---|---|---|
| `rt-circ-origin-any` | 1→1 | no hits | Bug |
| `rt-circ-origin-start` | no hits | same | Agrees |
| `rt-circ-origin-end` | 1→1 | no hits | Bug |
| `rt-circ-origin-within` | no hits | same | Agrees |
| `rt-circ-origin-equal` | no hits | same | Agrees |
| `rt-circ-origin-countOverlaps` | [1] | [0] | Bug |
| `rt-circ-origin-overlapsAny` | [true] | [false] | Bug |
| `rt-circ-origin-subsetByOverlaps` | 9/4; names=[origin]; label=[origin]; strand=[+] | empty; names=[]; label=[]; levels=[chrM] | Bug |
| `rt-circ-origin-nearest` | [1] | same | Agrees |
| `rt-circ-origin-precede` | [NA] | same | Agrees |
| `rt-circ-origin-follow` | [1] | same | Agrees |
| `rt-circ-origin-distance` | [6] | same | Agrees |
| `rt-circ-nearest-choice` | [2] | same | Agrees |
| `rt-circ-nearest-all` | 1→1, 1→2 | 1→2 | Bug |
| `rt-circ-nearest-origin-adjacent` | [1] | [2] | Bug |
| `rt-circ-origin-invert` | empty; names=[]; label=[]; levels=[chrM] | 9/4; names=[origin]; label=[origin]; strand=[+] | Bug |
| `rt-circ-origin-minoverlap` | 1→1 | no hits | Bug |
| `rt-circ-origin-opposite` | no hits | same | Agrees |
| `rt-circ-origin-ignore-strand` | 1→1 | no hits | Bug |
| `rt-circ-origin-linear-control` | no hits | same | Agrees |
| `rt-circ-origin-unknown-control` | no hits | same | Agrees |
| `rt-circ-origin-unknown-length` | no hits | same | Agrees |
| `rt-circ-shift-equal` | 1→1 | no hits | Bug |
| `rt-circ-shift-start` | 1→1 | no hits | Bug |
| `rt-circ-shift-within` | 1→1 | no hits | Bug |
| `rt-circ-zero-interior` | 1→1 | no hits | Bug |
| `rt-circ-zero-boundary` | no hits | same | Agrees |
| `rt-circ-multi-turn` | 1→1 | no hits | Bug |
| `rt-circ-self-drops` | 1→2, 2→1 | no hits | Bug |
| `rt-circ-self-redundant` | 1→2 | no hits | Bug |
| `rt-circ-origin-coverage` | chrM: runs=[1,0,1]; lengths=[2,6,2] | same | Agrees |
| `rt-circ-origin-trim` | 9/4; names=[origin]; label=[origin]; strand=[+] | same | Agrees |
| `rt-circ-origin-reduce` | 9/4; revmap=[[1]]; strand=[+] | same | Agrees |
