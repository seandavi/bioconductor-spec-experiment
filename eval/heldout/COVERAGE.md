# Heldout coverage

3510 cases. Cells = cases carrying the feature (required: >= 5 per family for the first seven).

| family | n | empty | zero_width | names | mcols | strand | seqinfo | na_arg | seqinfo_unused | seqinfo_unsorted | seqinfo_lengths | seqinfo_circular | seqinfo_genome | c_merge_seqinfo | trim_seqlengths | order_cross_seqlevel_strand | revmap | invalid | extreme |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| construct | 240 | 13 | 56 | 36 | 56 | 39 | 136 | 39 | 0 | 71 | 59 | 31 | 53 | 0 | 0 | 0 | 0 | 29 | 2 |
| accessors | 480 | 31 | 179 | 82 | 80 | 128 | 192 | 8 | 126 | 122 | 62 | 42 | 85 | 0 | 0 | 0 | 0 | 0 | 11 |
| metadata | 120 | 10 | 36 | 26 | 70 | 26 | 44 | 15 | 27 | 31 | 11 | 12 | 19 | 0 | 0 | 0 | 0 | 33 | 2 |
| seqinfo | 570 | 19 | 119 | 66 | 59 | 156 | 570 | 83 | 215 | 352 | 202 | 131 | 217 | 30 | 30 | 0 | 0 | 38 | 5 |
| subset_concat | 120 | 16 | 59 | 51 | 41 | 29 | 53 | 5 | 42 | 36 | 26 | 18 | 22 | 21 | 0 | 0 | 0 | 3 | 3 |
| intra | 300 | 19 | 99 | 41 | 54 | 81 | 104 | 29 | 69 | 67 | 36 | 25 | 53 | 0 | 0 | 0 | 0 | 57 | 3 |
| inter | 300 | 23 | 108 | 54 | 69 | 100 | 115 | 21 | 88 | 73 | 33 | 21 | 53 | 0 | 0 | 0 | 116 | 16 | 9 |
| setops | 360 | 26 | 178 | 118 | 104 | 146 | 174 | 12 | 107 | 122 | 61 | 63 | 82 | 0 | 0 | 0 | 0 | 15 | 7 |
| overlaps | 300 | 28 | 147 | 92 | 99 | 121 | 136 | 23 | 79 | 95 | 49 | 40 | 64 | 0 | 0 | 0 | 0 | 39 | 7 |
| nearest | 300 | 21 | 138 | 89 | 96 | 126 | 138 | 11 | 67 | 103 | 53 | 38 | 69 | 0 | 0 | 0 | 0 | 38 | 6 |
| coverage | 120 | 7 | 40 | 17 | 22 | 25 | 45 | 11 | 32 | 29 | 11 | 10 | 14 | 0 | 0 | 0 | 0 | 14 | 1 |
| ordering | 300 | 14 | 104 | 50 | 75 | 128 | 150 | 9 | 111 | 101 | 43 | 32 | 61 | 0 | 0 | 62 | 0 | 10 | 7 |

Legend: empty = a range input of length 0; zero_width = a width-0 range; names; mcols; strand = mixed
`+ - *` or `ignore.strand` given; seqinfo = multi-level / lengths / circular / genome; na_arg = NA in an
atomic argument; invalid = a deliberately invalid argument value was injected; extreme = coordinates near int32 limits.

## Cases per (op, class)

| op | class | n |
|---|---|---|
| GRanges | GRanges | 90 |
| IRanges | IRanges | 90 |
| Seqinfo | Seqinfo | 60 |
| [ | GRanges | 30 |
| [ | IRanges | 30 |
| c | GRanges | 30 |
| c | IRanges | 30 |
| countOverlaps | GRanges | 30 |
| countOverlaps | IRanges | 30 |
| coverage | GRanges | 60 |
| coverage | IRanges | 60 |
| disjoin | GRanges | 30 |
| disjoin | IRanges | 30 |
| distance | GRanges | 30 |
| distance | IRanges | 30 |
| duplicated | GRanges | 30 |
| duplicated | IRanges | 30 |
| end | GRanges | 30 |
| end | IRanges | 30 |
| findOverlaps | GRanges | 60 |
| findOverlaps | IRanges | 60 |
| flank | GRanges | 30 |
| flank | IRanges | 30 |
| follow | GRanges | 30 |
| follow | IRanges | 30 |
| gaps | GRanges | 30 |
| gaps | IRanges | 30 |
| genome | GRanges | 30 |
| genome | Seqinfo | 30 |
| genome<- | GRanges | 30 |
| genome<- | Seqinfo | 30 |
| intersect | GRanges | 30 |
| intersect | IRanges | 30 |
| isCircular | GRanges | 30 |
| isCircular | Seqinfo | 30 |
| isCircular<- | GRanges | 30 |
| isCircular<- | Seqinfo | 30 |
| length | GRanges | 30 |
| length | IRanges | 30 |
| mcols | GRanges | 30 |
| mcols | IRanges | 30 |
| mcols<- | GRanges | 30 |
| mcols<- | IRanges | 30 |
| merge | Seqinfo | 30 |
| names | GRanges | 30 |
| names | IRanges | 30 |
| names<- | GRanges | 30 |
| names<- | IRanges | 30 |
| narrow | GRanges | 30 |
| narrow | IRanges | 30 |
| nearest | GRanges | 60 |
| nearest | IRanges | 60 |
| order | GRanges | 30 |
| order | IRanges | 30 |
| overlapsAny | GRanges | 30 |
| overlapsAny | IRanges | 30 |
| pintersect | GRanges | 30 |
| pintersect | IRanges | 30 |
| precede | GRanges | 30 |
| precede | IRanges | 30 |
| promoters | GRanges | 30 |
| promoters | IRanges | 30 |
| psetdiff | GRanges | 30 |
| psetdiff | IRanges | 30 |
| punion | GRanges | 30 |
| punion | IRanges | 30 |
| range | GRanges | 30 |
| range | IRanges | 30 |
| ranges | GRanges | 30 |
| ranges | IRanges | 30 |
| rank | GRanges | 30 |
| rank | IRanges | 30 |
| reduce | GRanges | 60 |
| reduce | IRanges | 60 |
| resize | GRanges | 30 |
| resize | IRanges | 30 |
| seqinfo | GRanges | 30 |
| seqinfo<- | GRanges | 30 |
| seqlengths | GRanges | 30 |
| seqlengths | Seqinfo | 30 |
| seqlengths<- | GRanges | 30 |
| seqlengths<- | Seqinfo | 30 |
| seqlevels | GRanges | 30 |
| seqlevels | Seqinfo | 30 |
| seqlevels<- | GRanges | 30 |
| seqnames | GRanges | 30 |
| setdiff | GRanges | 30 |
| setdiff | IRanges | 30 |
| shift | GRanges | 30 |
| shift | IRanges | 30 |
| sort | GRanges | 30 |
| sort | IRanges | 30 |
| start | GRanges | 30 |
| start | IRanges | 30 |
| strand | GRanges | 30 |
| subsetByOverlaps | GRanges | 30 |
| subsetByOverlaps | IRanges | 30 |
| trim | GRanges | 30 |
| union | GRanges | 30 |
| union | IRanges | 30 |
| unique | GRanges | 30 |
| unique | IRanges | 30 |
| width | GRanges | 30 |
| width | IRanges | 30 |

Invalid-injected: 292 (8.3 %); extreme-coordinate: 63 (1.8 %).
