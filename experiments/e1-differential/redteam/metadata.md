# Red team: metadata

40 cases: 40 agree, 0 strict expected failures. Accepted divergences are pinned to an ADR.

Legend: coordinates are start/width; hits/indices are R's 1-based values; NA is missing.
Tables abbreviate outputs; the comparator checks complete metadata, strand, and seqinfo.
'Agrees' includes matching errors. 'ADR judgment' is unresolved, not permission to diverge.

| Case id | R | Python | View |
|---|---|---|---|
| `rt-meta-ir-construct-mcols-scalar` | 1/2,3/2,5/2; names=[r0,r1,r2]; score=[7,7,7] | same | Agrees |
| `rt-meta-ir-set-mcols-scalar` | 1/2,3/2,5/2; names=[r0,r1,r2]; score=[7,7,7] | same | Agrees |
| `rt-meta-ir-construct-mcols-divisor` | 1/2,3/2,5/2,7/2; names=[r0,r1,r2,r3]; score=[7,8,7,8] | same | Agrees |
| `rt-meta-ir-set-mcols-divisor` | 1/2,3/2,5/2,7/2; names=[r0,r1,r2,r3]; score=[7,8,7,8] | same | Agrees |
| `rt-meta-ir-construct-mcols-nondivisor` | 1/2,3/2,5/2; names=[r0,r1,r2]; score=[7,8,7] | same | Agrees |
| `rt-meta-ir-set-mcols-nondivisor` | 1/2,3/2,5/2; names=[r0,r1,r2]; score=[7,8,7] | same | Agrees |
| `rt-meta-ir-construct-mcols-zero` | error: trying to set a zero length metadata column on a non-zero length object | also errors | Agrees |
| `rt-meta-ir-set-mcols-zero` | error: trying to set a zero length metadata column on a non-zero length object | also errors | Agrees |
| `rt-meta-ir-flank-na-logical` | error: 'start' must be logical without NA's | also errors | Agrees |
| `rt-meta-ir-coverage-weight-column` | runs=[1,2,3]; lengths=[2,2,2] | same | Agrees |
| `rt-meta-gr-construct-mcols-scalar` | 1/2,3/2,5/2; names=[r0,r1,r2]; score=[7,7,7]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-set-mcols-scalar` | 1/2,3/2,5/2; names=[r0,r1,r2]; score=[7,7,7]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-construct-mcols-divisor` | 1/2,3/2,5/2,7/2; names=[r0,r1,r2,r3]; score=[7,8,7,8]; strand=[+,+,+,+] | same | Agrees |
| `rt-meta-gr-set-mcols-divisor` | 1/2,3/2,5/2,7/2; names=[r0,r1,r2,r3]; score=[7,8,7,8]; strand=[+,+,+,+] | same | Agrees |
| `rt-meta-gr-construct-mcols-nondivisor` | 1/2,3/2,5/2; names=[r0,r1,r2]; score=[7,8,7]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-set-mcols-nondivisor` | 1/2,3/2,5/2; names=[r0,r1,r2]; score=[7,8,7]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-construct-mcols-zero` | error: trying to set a zero length metadata column on a non-zero length object | also errors | Agrees |
| `rt-meta-gr-set-mcols-zero` | error: trying to set a zero length metadata column on a non-zero length object | also errors | Agrees |
| `rt-meta-gr-flank-na-logical` | error: 'start' must be logical without NA's | also errors | Agrees |
| `rt-meta-gr-coverage-weight-column` | chr1: runs=[1,2,3]; lengths=[2,2,2] | same | Agrees |
| `rt-meta-gr-reserved-start` | 1/2,3/2,5/2; names=[r0,r1,r2]; start=[1,2,3]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-reserved-seqnames` | 1/2,3/2,5/2; names=[r0,r1,r2]; seqnames=[1,2,3]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-reserved-strand` | 1/2,3/2,5/2; names=[r0,r1,r2]; strand=[1,2,3]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-reserved-width` | 1/2,3/2,5/2; names=[r0,r1,r2]; width=[1,2,3]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-nonreserved-control` | 1/2,3/2,5/2; names=[r0,r1,r2]; hit=[false,false,false]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-pintersect-existing-hit-keep` | 2/1,3/0,5/2; names=[r0,r1,r2]; hit=[true,false,true]; label=[a,b,c]; strand=[+,+,+] | same | Agrees |
| `rt-meta-gr-pintersect-existing-hit-drop` | 2/1,5/2; names=[r0,r2]; hit=[false,false]; label=[a,c]; strand=[+,+] | same | Agrees |
| `rt-meta-gr-coverage-shift-column` | chr1: runs=[1,0,1,0,1]; lengths=[2,1,2,1,2] | same | Agrees |
| `rt-meta-gr-coverage-na-column` | chr1: runs=[1,NA,3]; lengths=[2,2,2] | same | Agrees |
| `rt-meta-gr-strand-text-na` | error: strand values must be in '+' '-' '*' | also errors | Agrees |
| `rt-meta-gr-strand-numeric` | error: strand values must be in '1' '-1' 'NA' | also errors | Agrees |
| `rt-meta-gr-strand-typed-missing` | 1/2,3/2,5/2; names=[r0,r1,r2]; strand=[+,*,-] | same | Agrees |
| `rt-meta-seqinfo-numeric-circular` | error: bad supplied 'isCircular' vector | also errors | Agrees |
| `rt-meta-seqinfo-text-circular` | error: bad supplied 'isCircular' vector | also errors | Agrees |
| `rt-meta-seqinfo-empty` | levels=[]; circular=[] | same | Agrees |
| `rt-meta-seqinfo-numeric-names` | error: supplied 'seqlevels' must be a character vector | also errors | Agrees |
| `rt-meta-gr-set-seqinfo-reorder` | error: The seqlevels in the supplied Seqinfo object are not in the same order as the | also errors | Agrees |
| `rt-meta-gr-set-seqinfo-drop-unused` | error: The seqlevels in the supplied Seqinfo object are not the same as the | also errors | Agrees |
| `rt-meta-gr-set-seqinfo-drop-used` | error: The seqlevels in the supplied Seqinfo object are not the same as the | also errors | Agrees |
| `rt-meta-gr-recycle-range-names` | 1/2,1/2; names=[only,only]; strand=[*,*]; levels=[a,b] | same | Agrees |
