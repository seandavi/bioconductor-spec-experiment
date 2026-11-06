# Fidelity to R: evaluation report

R: `R 4.6.0; IRanges 2.47.2; GenomicRanges 1.65.1; S4Vectors 0.51.6; Seqinfo 1.3.0; BiocGenerics 0.59.12; jsonlite 2.0.0`
Oracle: 9650 cases, 9650 scored, 0 excluded. Full R run: 9577 distinct cases in 21.9 s (4 workers).

## Pool and reconciliation

| suite | read/generated | converted | unconvertible | R-excluded | scored |
|---|---|---|---|---|---|
| e1 | 4515 | 4504 | 11 | 0 | 4504 |
| e2 | 921 | 921 | 0 | 0 | 921 |
| e4 | 715 | 715 | 0 | 0 | 715 |
| heldout | 3510 | 3510 | 0 | 0 | 3510 |

read = converted + unconvertible; converted = R-excluded + scored (asserted when this report is built).

Unconvertible (listed in `cases/unconvertible.jsonl`):

| suite | reason | n |
|---|---|---|
| e1 | literal:source-multi-key | 7 |
| e1 | op-out-of-vocab:as.vector | 1 |
| e1 | op-out-of-vocab:runLength | 1 |
| e1 | op-out-of-vocab:runValue | 2 |

## Headline: pass rate per port x source

Pass rate = pass / scored; unsupported, adapter errors and port errors stay in the denominator. **heldout is the headline** (no port has seen it). `*` = the port's home suite (it was built against it).

| port | heldout | e1 | e1-redteam | e2 | e4 | all |
|---|---|---|---|---|---|---|
| e0 | 2781/3510 (79.2%) [mism155 unsup542 perr32] | 3227/4385 (73.6%) [mism935 unsup96 perr127] | 69/119 (58.0%) [mism31 unsup1 perr18] | 836/921 (90.8%) [mism35 unsup33 perr17] | 674/715 (94.3%) [mism30 unsup5 perr6] | 7587/9650 (78.6%) |
| e0b | 2510/3510 (71.5%) [mism178 unsup743 perr79] | 2618/4385 (59.7%) [mism800 unsup505 perr462] | 60/119 (50.4%) [mism25 unsup15 perr19] | 765/921 (83.1%) [mism59 unsup51 perr46] | 603/715 (84.3%) [mism68 unsup12 perr32] | 6556/9650 (67.9%) |
| e1 | 2987/3510 (85.1%) [mism77 unsup436 perr10] | *3940/4385 (89.9%) [mism370 unsup73 perr2] | *110/119 (92.4%) [mism9] | 876/921 (95.1%) [mism19 unsup26] | 706/715 (98.7%) [mism6 perr3] | 8619/9650 (89.3%) |
| e2 | 2896/3510 (82.5%) [mism117 unsup478 perr19] | 3132/4385 (71.4%) [mism649 unsup235 perr369] | 84/119 (70.6%) [mism27 unsup2 perr6] | *891/921 (96.7%) [mism15 unsup15] | 685/715 (95.8%) [mism23 unsup4 perr3] | 7688/9650 (79.7%) |
| e3 | 2573/3510 (73.3%) [mism144 unsup739 perr54] | 2918/4385 (66.5%) [mism532 unsup698 perr237] | 61/119 (51.3%) [mism25 unsup17 perr16] | 804/921 (87.3%) [mism32 unsup64 perr21] | 639/715 (89.4%) [mism39 unsup21 perr16] | 6995/9650 (72.5%) |
| e4 | 2645/3510 (75.4%) [mism123 unsup709 perr33] | 2954/4385 (67.4%) [mism885 unsup378 perr168] | 61/119 (51.3%) [mism30 unsup12 perr16] | 814/921 (88.4%) [mism31 unsup62 perr14] | *698/715 (97.6%) [mism8 unsup1 perr8] | 7172/9650 (74.3%) |
| biocpy | 2160/3510 (61.5%) [mism423 unsup739 perr188] | 1929/4385 (44.0%) [mism1364 unsup517 perr575] | 41/119 (34.5%) [mism38 unsup10 perr30] | 564/921 (61.2%) [mism181 unsup121 perr55] | 414/715 (57.9%) [mism168 unsup71 perr62] | 5108/9650 (52.9%) |

Cell: `pass/scored (rate) [mism = mismatch, adapt = adapter_error, unsup = unsupported, perr = port_error_unexpected]`.

## e0: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 226/240 (94.2%) [mism9 perr5] | 311/399 (77.9%) [mism19 perr69] | 3/10 (30.0%) [mism4 unsup1 perr2] | 50/60 (83.3%) [mism3 unsup2 perr5] | 41/46 (89.1%) [mism5] |
| accessors | 407/480 (84.8%) [mism6 unsup67] | 32/36 (88.9%) [mism3 unsup1] | - | 24/27 (88.9%) [mism3] | 14/15 (93.3%) [unsup1] |
| metadata | 104/120 (86.7%) [mism11 perr5] | 21/24 (87.5%) [mism2 unsup1] | 9/21 (42.9%) [perr12] | 9/11 (81.8%) [perr2] | 2/2 (100.0%) |
| seqinfo | 315/570 (55.3%) [unsup255] | 107/114 (93.9%) [unsup7] | 6/8 (75.0%) [mism2] | 30/50 (60.0%) [mism1 unsup19] | 25/25 (100.0%) |
| subset_concat | 78/120 (65.0%) [mism1 unsup40 perr1] | 351/462 (76.0%) [mism33 unsup78] | 1/1 (100.0%) | 50/60 (83.3%) [unsup10] | 36/37 (97.3%) [unsup1] |
| intra | 283/300 (94.3%) [mism4 perr13] | 321/363 (88.4%) [mism14 perr28] | 2/2 (100.0%) | 124/126 (98.4%) [perr2] | 82/88 (93.2%) [mism5 perr1] |
| inter | 240/300 (80.0%) [mism14 unsup39 perr7] | 340/405 (84.0%) [mism61 perr4] | 12/12 (100.0%) | 91/101 (90.1%) [mism6 perr4] | 101/104 (97.1%) [mism3] |
| setops | 338/360 (93.9%) [mism16 unsup6] | 397/441 (90.0%) [mism36 perr8] | 4/4 (100.0%) | 102/107 (95.3%) [mism5] | 106/108 (98.1%) [mism1 perr1] |
| overlaps | 273/300 (91.0%) [mism27] | 539/848 (63.6%) [mism304 unsup3 perr2] | 8/23 (34.8%) [mism15] | 185/195 (94.9%) [mism10] | 109/113 (96.5%) [mism3 unsup1] |
| nearest | 256/300 (85.3%) [mism43 perr1] | 528/582 (90.7%) [mism49 perr5] | 10/15 (66.7%) [mism5] | 72/76 (94.7%) [mism3 perr1] | 62/67 (92.5%) [mism4 perr1] |
| coverage | 70/120 (58.3%) [mism10 unsup40] | 117/511 (22.9%) [mism379 unsup6 perr9] | 3/12 (25.0%) [mism5 perr4] | 44/51 (86.3%) [mism3 unsup1 perr3] | 41/54 (75.9%) [mism8 unsup2 perr3] |
| ordering | 191/300 (63.7%) [mism14 unsup95] | 163/200 (81.5%) [mism35 perr2] | 11/11 (100.0%) | 55/57 (96.5%) [mism1 unsup1] | 55/56 (98.2%) [mism1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 256, 'exact': 693, 'R2': 197, 'R3': 40}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 5 | 24 | 22 | 10 | 0 | 63 | 13 | 20 |
| e1 | 1 | 268 | 19 | 26 | 1 | 169 | 249 | 439 |
| e1-redteam | 0 | 19 | 4 | 1 | 0 | 6 | 2 | 4 |
| e2 | 3 | 2 | 10 | 0 | 1 | 8 | 3 | 8 |
| e4 | 2 | 10 | 2 | 0 | 0 | 10 | 1 | 5 |

e0 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 321/443 (72.5%) [mism33 unsup84 perr5] |
| empty | 185/223 (83.0%) [mism14 unsup19 perr5] |
| mcols | 627/804 (78.0%) [mism48 unsup126 perr3] |
| multi_seqlevel | 1203/1645 (73.1%) [mism106 unsup323 perr13] |
| na_arg | 54/222 (24.3%) [mism56 unsup96 perr16] |
| names | 538/694 (77.5%) [mism49 unsup102 perr5] |
| r_error | 292/440 (66.4%) [mism63 unsup85] |
| seqlengths | 407/593 (68.6%) [mism45 unsup140 perr1] |
| strand | 884/1187 (74.5%) [mism82 unsup211 perr10] |
| zero_width | 969/1210 (80.1%) [mism62 unsup171 perr8] |

## e0b: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 219/240 (91.2%) [mism9 unsup7 perr5] | 273/399 (68.4%) [mism15 unsup31 perr80] | 6/10 (60.0%) [mism1 unsup1 perr2] | 48/60 (80.0%) [mism2 unsup5 perr5] | 40/46 (87.0%) [mism3 unsup2 perr1] |
| accessors | 402/480 (83.8%) [mism11 unsup67] | 31/36 (86.1%) [mism5] | - | 23/27 (85.2%) [mism4] | 15/15 (100.0%) |
| metadata | 115/120 (95.8%) [perr5] | 24/24 (100.0%) | 9/21 (42.9%) [perr12] | 9/11 (81.8%) [perr2] | 2/2 (100.0%) |
| seqinfo | 315/570 (55.3%) [unsup255] | 107/114 (93.9%) [unsup7] | 6/8 (75.0%) [mism2] | 30/50 (60.0%) [mism1 unsup19] | 25/25 (100.0%) |
| subset_concat | 57/120 (47.5%) [unsup43 perr20] | 225/462 (48.7%) [mism1 unsup78 perr158] | 0/1 (0.0%) [mism1] | 37/60 (61.7%) [unsup11 perr12] | 33/37 (89.2%) [mism1 perr3] |
| intra | 228/300 (76.0%) [mism6 unsup63 perr3] | 301/363 (82.9%) [mism15 unsup19 perr28] | 0/2 (0.0%) [unsup2] | 122/126 (96.8%) [unsup2 perr2] | 79/88 (89.8%) [mism6 unsup2 perr1] |
| inter | 185/300 (61.7%) [mism33 unsup79 perr3] | 196/405 (48.4%) [mism122 unsup80 perr7] | 9/12 (75.0%) [unsup3] | 80/101 (79.2%) [mism13 unsup5 perr3] | 87/104 (83.7%) [mism12 unsup5] |
| setops | 282/360 (78.3%) [mism44 unsup13 perr21] | 321/441 (72.8%) [mism34 unsup10 perr76] | 2/4 (50.0%) [unsup2] | 78/107 (72.9%) [mism18 unsup1 perr10] | 82/108 (75.9%) [mism19 perr7] |
| overlaps | 275/300 (91.7%) [mism15 unsup10] | 501/848 (59.1%) [mism287 unsup58 perr2] | 8/23 (34.8%) [mism13 unsup2] | 187/195 (95.9%) [mism8] | 108/113 (95.6%) [mism4 unsup1] |
| nearest | 231/300 (77.0%) [mism52 perr17] | 396/582 (68.0%) [mism108 perr78] | 10/15 (66.7%) [mism5] | 60/76 (78.9%) [mism10 perr6] | 43/67 (64.2%) [mism11 perr13] |
| coverage | 73/120 (60.8%) [mism5 unsup42] | 111/511 (21.7%) [mism195 unsup180 perr25] | 3/12 (25.0%) [mism3 unsup2 perr4] | 44/51 (86.3%) [mism2 perr5] | 38/54 (70.4%) [mism11 perr5] |
| ordering | 128/300 (42.7%) [mism3 unsup164 perr5] | 132/200 (66.0%) [mism18 unsup42 perr8] | 7/11 (63.6%) [unsup3 perr1] | 47/57 (82.5%) [mism1 unsup8 perr1] | 51/56 (91.1%) [mism1 unsup2 perr2] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 183, 'exact': 680, 'R2': 221, 'R3': 46}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 7 | 65 | 19 | 12 | 7 | 43 | 0 | 30 |
| e1 | 2 | 274 | 15 | 30 | 58 | 117 | 1 | 324 |
| e1-redteam | 0 | 15 | 3 | 1 | 0 | 5 | 0 | 4 |
| e2 | 3 | 22 | 9 | 0 | 5 | 6 | 1 | 14 |
| e4 | 2 | 36 | 1 | 0 | 5 | 12 | 1 | 15 |

e0b heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 289/443 (65.2%) [mism32 unsup110 perr12] |
| empty | 170/223 (76.2%) [mism13 unsup31 perr9] |
| mcols | 543/804 (67.5%) [mism58 unsup171 perr32] |
| multi_seqlevel | 1076/1645 (65.4%) [mism103 unsup426 perr40] |
| na_arg | 81/222 (36.5%) [mism29 unsup109 perr3] |
| names | 480/694 (69.2%) [mism55 unsup133 perr26] |
| r_error | 288/440 (65.5%) [mism43 unsup109] |
| seqlengths | 376/593 (63.4%) [mism34 unsup170 perr13] |
| strand | 787/1187 (66.3%) [mism84 unsup285 perr31] |
| zero_width | 847/1210 (70.0%) [mism83 unsup241 perr39] |

## e1: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 227/240 (94.6%) [mism10 perr3] | 399/399 (100.0%) | 10/10 (100.0%) | 57/60 (95.0%) [mism3] | 46/46 (100.0%) |
| accessors | 407/480 (84.8%) [mism6 unsup67] | 32/36 (88.9%) [mism4] | - | 24/27 (88.9%) [mism3] | 15/15 (100.0%) |
| metadata | 120/120 (100.0%) | 24/24 (100.0%) | 21/21 (100.0%) | 11/11 (100.0%) | 2/2 (100.0%) |
| seqinfo | 416/570 (73.0%) [mism1 unsup149 perr4] | 114/114 (100.0%) | 8/8 (100.0%) | 35/50 (70.0%) [unsup15] | 25/25 (100.0%) |
| subset_concat | 78/120 (65.0%) [mism1 unsup40 perr1] | 388/462 (84.0%) [mism1 unsup73] | 1/1 (100.0%) | 50/60 (83.3%) [unsup10] | 37/37 (100.0%) |
| intra | 298/300 (99.3%) [mism1 perr1] | 362/363 (99.7%) [perr1] | 2/2 (100.0%) | 126/126 (100.0%) | 85/88 (96.6%) [perr3] |
| inter | 261/300 (87.0%) [unsup39] | 405/405 (100.0%) | 12/12 (100.0%) | 101/101 (100.0%) | 102/104 (98.1%) [mism2] |
| setops | 336/360 (93.3%) [mism18 unsup6] | 441/441 (100.0%) | 4/4 (100.0%) | 102/107 (95.3%) [mism5] | 107/108 (99.1%) [mism1] |
| overlaps | 294/300 (98.0%) [mism6] | 835/848 (98.5%) [mism13] | 22/23 (95.7%) [mism1] | 188/195 (96.4%) [mism7] | 111/113 (98.2%) [mism2] |
| nearest | 272/300 (90.7%) [mism28] | 582/582 (100.0%) | 12/15 (80.0%) [mism3] | 76/76 (100.0%) | 67/67 (100.0%) |
| coverage | 76/120 (63.3%) [mism3 unsup40 perr1] | 159/511 (31.1%) [mism351 perr1] | 7/12 (58.3%) [mism5] | 51/51 (100.0%) | 54/54 (100.0%) |
| ordering | 202/300 (67.3%) [mism3 unsup95] | 199/200 (99.5%) [mism1] | 11/11 (100.0%) | 55/57 (96.5%) [mism1 unsup1] | 55/56 (98.2%) [mism1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 11, 'R2': 14, 'exact': 446, 'R3': 10}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 2 | 23 | 19 | 10 | 1 | 3 | 3 | 18 |
| e1 | 0 | 3 | 13 | 0 | 0 | 3 | 180 | 341 |
| e1-redteam | 0 | 2 | 1 | 0 | 0 | 1 | 3 | 5 |
| e2 | 1 | 1 | 9 | 0 | 2 | 2 | 0 | 4 |
| e4 | 0 | 2 | 1 | 0 | 0 | 2 | 0 | 1 |

e1 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 368/443 (83.1%) [mism17 unsup57 perr1] |
| empty | 201/223 (90.1%) [mism3 unsup15 perr4] |
| mcols | 671/804 (83.5%) [mism24 unsup108 perr1] |
| multi_seqlevel | 1358/1645 (82.6%) [mism50 unsup231 perr6] |
| na_arg | 158/222 (71.2%) [mism3 unsup59 perr2] |
| names | 575/694 (82.9%) [mism38 unsup79 perr2] |
| r_error | 363/440 (82.5%) [mism3 unsup74] |
| seqlengths | 471/593 (79.4%) [mism20 unsup101 perr1] |
| strand | 1001/1187 (84.3%) [mism43 unsup140 perr3] |
| zero_width | 1042/1210 (86.1%) [mism32 unsup134 perr2] |

## e2: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 235/240 (97.9%) [perr5] | 337/399 (84.5%) [mism27 perr35] | 6/10 (60.0%) [mism3 perr1] | 59/60 (98.3%) [unsup1] | 44/46 (95.7%) [mism1 perr1] |
| accessors | 403/480 (84.0%) [mism10 unsup67] | 31/36 (86.1%) [mism5] | - | 24/27 (88.9%) [mism3] | 15/15 (100.0%) |
| metadata | 120/120 (100.0%) | 24/24 (100.0%) | 21/21 (100.0%) | 11/11 (100.0%) | 2/2 (100.0%) |
| seqinfo | 449/570 (78.8%) [mism5 unsup113 perr3] | 113/114 (99.1%) [mism1] | 7/8 (87.5%) [mism1] | 46/50 (92.0%) [unsup4] | 25/25 (100.0%) |
| subset_concat | 75/120 (62.5%) [mism1 unsup43 perr1] | 351/462 (76.0%) [mism33 unsup78] | 0/1 (0.0%) [mism1] | 51/60 (85.0%) [unsup9] | 36/37 (97.3%) [mism1] |
| intra | 239/300 (79.7%) [mism10 unsup44 perr7] | 308/363 (84.8%) [mism33 unsup14 perr8] | 2/2 (100.0%) | 126/126 (100.0%) | 85/88 (96.6%) [mism2 perr1] |
| inter | 242/300 (80.7%) [mism11 unsup47] | 303/405 (74.8%) [mism60 unsup42] | 12/12 (100.0%) | 101/101 (100.0%) | 102/104 (98.1%) [mism2] |
| setops | 334/360 (92.8%) [mism13 unsup13] | 388/441 (88.0%) [mism44 unsup1 perr8] | 3/4 (75.0%) [mism1] | 107/107 (100.0%) | 106/108 (98.1%) [mism1 perr1] |
| overlaps | 271/300 (90.3%) [mism13 unsup16] | 494/848 (58.3%) [mism252 unsup100 perr2] | 8/23 (34.8%) [mism13 unsup2] | 188/195 (96.4%) [mism7] | 107/113 (94.7%) [mism2 unsup4] |
| nearest | 256/300 (85.3%) [mism44] | 508/582 (87.3%) [mism49 perr25] | 10/15 (66.7%) [mism5] | 73/76 (96.1%) [mism3] | 63/67 (94.0%) [mism4] |
| coverage | 75/120 (62.5%) [mism2 unsup40 perr3] | 128/511 (25.0%) [mism94 perr289] | 4/12 (33.3%) [mism3 perr5] | 50/51 (98.0%) [mism1] | 45/54 (83.3%) [mism9] |
| ordering | 197/300 (65.7%) [mism8 unsup95] | 147/200 (73.5%) [mism51 perr2] | 11/11 (100.0%) | 55/57 (96.5%) [mism1 unsup1] | 55/56 (98.2%) [mism1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 293, 'exact': 338, 'R2': 156, 'R3': 44}

| source | kind | length | names | not_in_universe | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|
| heldout | 8 | 7 | 19 | 11 | 57 | 1 | 14 |
| e1 | 2 | 191 | 15 | 29 | 222 | 1 | 190 |
| e1-redteam | 0 | 17 | 3 | 1 | 6 | 1 | 3 |
| e2 | 3 | 0 | 9 | 0 | 2 | 1 | 0 |
| e4 | 2 | 9 | 1 | 0 | 6 | 1 | 4 |

e2 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 361/443 (81.5%) [mism22 unsup58 perr2] |
| empty | 190/223 (85.2%) [mism9 unsup21 perr3] |
| mcols | 644/804 (80.1%) [mism26 unsup133 perr1] |
| multi_seqlevel | 1325/1645 (80.5%) [mism74 unsup238 perr8] |
| na_arg | 120/222 (54.1%) [mism32 unsup60 perr10] |
| names | 552/694 (79.5%) [mism44 unsup96 perr2] |
| r_error | 311/440 (70.7%) [mism57 unsup72] |
| seqlengths | 455/593 (76.7%) [mism30 unsup107 perr1] |
| strand | 933/1187 (78.6%) [mism60 unsup187 perr7] |
| zero_width | 993/1210 (82.1%) [mism52 unsup164 perr1] |

## e3: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 214/240 (89.2%) [mism9 unsup6 perr11] | 290/399 (72.7%) [mism24 perr85] | 4/10 (40.0%) [mism3 unsup1 perr2] | 47/60 (78.3%) [mism2 unsup5 perr6] | 38/46 (82.6%) [mism3 perr5] |
| accessors | 399/480 (83.1%) [mism14 unsup67] | 28/36 (77.8%) [mism8] | - | 21/27 (77.8%) [mism5 unsup1] | 14/15 (93.3%) [mism1] |
| metadata | 115/120 (95.8%) [perr5] | 24/24 (100.0%) | 9/21 (42.9%) [perr12] | 9/11 (81.8%) [perr2] | 2/2 (100.0%) |
| seqinfo | 314/570 (55.1%) [unsup256] | 107/114 (93.9%) [unsup7] | 5/8 (62.5%) [unsup3] | 29/50 (58.0%) [unsup21] | 25/25 (100.0%) |
| subset_concat | 74/120 (61.7%) [mism1 unsup45] | 351/462 (76.0%) [mism33 unsup78] | 0/1 (0.0%) [mism1] | 49/60 (81.7%) [unsup11] | 36/37 (97.3%) [mism1] |
| intra | 213/300 (71.0%) [mism6 unsup72 perr9] | 272/363 (74.9%) [mism35 unsup21 perr35] | 0/2 (0.0%) [mism2] | 120/126 (95.2%) [unsup4 perr2] | 78/88 (88.6%) [mism5 unsup4 perr1] |
| inter | 232/300 (77.3%) [mism11 unsup50 perr7] | 313/405 (77.3%) [mism38 unsup43 perr11] | 12/12 (100.0%) | 96/101 (95.0%) [unsup1 perr4] | 101/104 (97.1%) [mism2 unsup1] |
| setops | 303/360 (84.2%) [mism41 unsup13 perr3] | 368/441 (83.4%) [mism54 unsup10 perr9] | 2/4 (50.0%) [unsup2] | 94/107 (87.9%) [mism11 unsup1 perr1] | 95/108 (88.0%) [mism12 perr1] |
| overlaps | 273/300 (91.0%) [mism11 unsup16] | 496/848 (58.5%) [mism250 unsup100 perr2] | 8/23 (34.8%) [mism13 unsup2] | 188/195 (96.4%) [mism7] | 107/113 (94.7%) [mism2 unsup4] |
| nearest | 239/300 (79.7%) [mism47 perr14] | 457/582 (78.5%) [mism49 perr76] | 10/15 (66.7%) [mism5] | 67/76 (88.2%) [mism6 perr3] | 56/67 (83.6%) [mism5 perr6] |
| coverage | 43/120 (35.8%) [unsup77] | 70/511 (13.7%) [mism14 unsup416 perr11] | 3/12 (25.0%) [mism1 unsup7 perr1] | 35/51 (68.6%) [unsup14 perr2] | 34/54 (63.0%) [mism7 unsup12 perr1] |
| ordering | 154/300 (51.3%) [mism4 unsup137 perr5] | 142/200 (71.0%) [mism27 unsup23 perr8] | 8/11 (72.7%) [unsup2 perr1] | 49/57 (86.0%) [mism1 unsup6 perr1] | 53/56 (94.6%) [mism1 perr2] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 227, 'exact': 339, 'R2': 160, 'R3': 46}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 7 | 41 | 22 | 12 | 0 | 46 | 1 | 16 |
| e1 | 2 | 197 | 18 | 30 | 0 | 161 | 1 | 123 |
| e1-redteam | 0 | 15 | 3 | 1 | 0 | 7 | 0 | 2 |
| e2 | 3 | 11 | 10 | 0 | 1 | 3 | 1 | 3 |
| e4 | 2 | 19 | 2 | 0 | 0 | 10 | 1 | 5 |

e3 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 296/443 (66.8%) [mism28 unsup110 perr9] |
| empty | 175/223 (78.5%) [mism11 unsup33 perr4] |
| mcols | 590/804 (73.4%) [mism39 unsup166 perr9] |
| multi_seqlevel | 1107/1645 (67.3%) [mism95 unsup423 perr20] |
| na_arg | 71/222 (32.0%) [mism29 unsup108 perr14] |
| names | 492/694 (70.9%) [mism56 unsup138 perr8] |
| r_error | 277/440 (63.0%) [mism46 unsup117] |
| seqlengths | 386/593 (65.1%) [mism33 unsup170 perr4] |
| strand | 808/1187 (68.1%) [mism75 unsup286 perr18] |
| zero_width | 883/1210 (73.0%) [mism71 unsup241 perr15] |

## e4: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 217/240 (90.4%) [mism9 unsup6 perr8] | 299/399 (74.9%) [mism19 perr81] | 4/10 (40.0%) [mism3 unsup1 perr2] | 49/60 (81.7%) [mism1 unsup5 perr5] | 45/46 (97.8%) [perr1] |
| accessors | 402/480 (83.8%) [mism11 unsup67] | 31/36 (86.1%) [mism5] | - | 23/27 (85.2%) [mism4] | 15/15 (100.0%) |
| metadata | 85/120 (70.8%) [unsup30 perr5] | 21/24 (87.5%) [unsup3] | 9/21 (42.9%) [perr12] | 8/11 (72.7%) [unsup1 perr2] | 2/2 (100.0%) |
| seqinfo | 314/570 (55.1%) [unsup256] | 107/114 (93.9%) [unsup7] | 5/8 (62.5%) [unsup3] | 29/50 (58.0%) [unsup21] | 25/25 (100.0%) |
| subset_concat | 75/120 (62.5%) [mism1 unsup43 perr1] | 351/462 (76.0%) [mism33 unsup78] | 1/1 (100.0%) | 49/60 (81.7%) [unsup11] | 37/37 (100.0%) |
| intra | 263/300 (87.7%) [mism4 unsup22 perr11] | 287/363 (79.1%) [mism27 unsup13 perr36] | 0/2 (0.0%) [mism2] | 124/126 (98.4%) [perr2] | 84/88 (95.5%) [perr4] |
| inter | 233/300 (77.7%) [mism19 unsup48] | 256/405 (63.2%) [mism75 unsup72 perr2] | 10/12 (83.3%) [mism2] | 90/101 (89.1%) [mism7 unsup4] | 102/104 (98.1%) [mism2] |
| setops | 317/360 (88.1%) [mism23 unsup13 perr7] | 369/441 (83.7%) [mism46 unsup18 perr8] | 2/4 (50.0%) [unsup2] | 97/107 (90.7%) [mism6 unsup1 perr3] | 107/108 (99.1%) [mism1] |
| overlaps | 277/300 (92.3%) [mism11 unsup12] | 501/848 (59.1%) [mism268 unsup77 perr2] | 8/23 (34.8%) [mism13 unsup2] | 188/195 (96.4%) [mism7] | 111/113 (98.2%) [mism2] |
| nearest | 235/300 (78.3%) [mism39 unsup25 perr1] | 447/582 (76.8%) [mism39 unsup70 perr26] | 8/15 (53.3%) [mism5 unsup2] | 58/76 (76.3%) [mism5 unsup12 perr1] | 63/67 (94.0%) [mism3 perr1] |
| coverage | 72/120 (60.0%) [mism5 unsup43] | 153/511 (29.9%) [mism347 perr11] | 5/12 (41.7%) [mism5 perr2] | 49/51 (96.1%) [mism1 perr1] | 52/54 (96.3%) [perr2] |
| ordering | 155/300 (51.7%) [mism1 unsup144] | 132/200 (66.0%) [mism26 unsup40 perr2] | 9/11 (81.8%) [unsup2] | 50/57 (87.7%) [unsup7] | 55/56 (98.2%) [unsup1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 192, 'exact': 673, 'R2': 174, 'R3': 38}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 8 | 35 | 19 | 8 | 0 | 37 | 1 | 16 |
| e1 | 2 | 453 | 15 | 26 | 8 | 147 | 1 | 244 |
| e1-redteam | 0 | 21 | 3 | 1 | 0 | 6 | 1 | 2 |
| e2 | 3 | 6 | 9 | 0 | 1 | 1 | 1 | 10 |
| e4 | 2 | 2 | 1 | 0 | 0 | 1 | 1 | 1 |

e4 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 300/443 (67.7%) [mism25 unsup114 perr4] |
| empty | 176/223 (78.9%) [mism11 unsup30 perr6] |
| mcols | 605/804 (75.2%) [mism34 unsup162 perr3] |
| multi_seqlevel | 1128/1645 (68.6%) [mism77 unsup432 perr8] |
| na_arg | 70/222 (31.5%) [mism36 unsup109 perr7] |
| names | 503/694 (72.5%) [mism51 unsup136 perr4] |
| r_error | 287/440 (65.2%) [mism37 unsup116] |
| seqlengths | 388/593 (65.4%) [mism28 unsup175 perr2] |
| strand | 817/1187 (68.8%) [mism64 unsup300 perr6] |
| zero_width | 918/1210 (75.9%) [mism56 unsup228 perr8] |

## biocpy: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 149/240 (62.1%) [mism17 unsup55 perr19] | 242/399 (60.7%) [mism24 unsup40 perr93] | 6/10 (60.0%) [mism2 perr2] | 36/60 (60.0%) [mism3 unsup13 perr8] | 27/46 (58.7%) [mism6 unsup7 perr6] |
| accessors | 400/480 (83.3%) [mism13 unsup67] | 28/36 (77.8%) [mism8] | - | 23/27 (85.2%) [mism4] | 13/15 (86.7%) [mism2] |
| metadata | 104/120 (86.7%) [mism6 perr10] | 20/24 (83.3%) [mism4] | 4/21 (19.0%) [perr17] | 9/11 (81.8%) [perr2] | 1/2 (50.0%) [mism1] |
| seqinfo | 467/570 (81.9%) [mism29 unsup62 perr12] | 37/114 (32.5%) [mism76 unsup1] | 5/8 (62.5%) [mism3] | 30/50 (60.0%) [mism13 unsup6 perr1] | 24/25 (96.0%) [mism1] |
| subset_concat | 53/120 (44.2%) [mism16 unsup43 perr8] | 303/462 (65.6%) [mism68 unsup78 perr13] | 0/1 (0.0%) [mism1] | 41/60 (68.3%) [mism8 unsup7 perr4] | 24/37 (64.9%) [mism11 perr2] |
| intra | 226/300 (75.3%) [mism7 unsup61 perr6] | 306/363 (84.3%) [mism35 unsup19 perr3] | 0/2 (0.0%) [mism2] | 118/126 (93.7%) [mism5 unsup2 perr1] | 73/88 (83.0%) [mism7 unsup2 perr6] |
| inter | 242/300 (80.7%) [mism17 unsup27 perr14] | 302/405 (74.6%) [mism90 perr13] | 8/12 (66.7%) [mism3 perr1] | 93/101 (92.1%) [mism6 perr2] | 93/104 (89.4%) [mism5 perr6] |
| setops | 113/360 (31.4%) [mism62 unsup157 perr28] | 108/441 (24.5%) [mism154 unsup163 perr16] | 1/4 (25.0%) [mism1 unsup2] | 45/107 (42.1%) [mism16 unsup41 perr5] | 36/108 (33.3%) [mism22 unsup36 perr14] |
| overlaps | 87/300 (29.0%) [mism98 unsup96 perr19] | 275/848 (32.4%) [mism407 unsup95 perr71] | 7/23 (30.4%) [mism10 unsup4 perr2] | 102/195 (52.3%) [mism53 unsup20 perr20] | 63/113 (55.8%) [mism34 unsup8 perr8] |
| nearest | 187/300 (62.3%) [mism49 unsup42 perr22] | 217/582 (37.3%) [mism215 unsup65 perr85] | 8/15 (53.3%) [mism5 unsup2] | 43/76 (56.6%) [mism16 unsup13 perr4] | 36/67 (53.7%) [mism17 unsup1 perr13] |
| coverage | 10/120 (8.3%) [mism79 unsup19 perr12] | 28/511 (5.5%) [mism204 perr279] | 0/12 (0.0%) [mism5 perr7] | 5/51 (9.8%) [mism39 perr7] | 3/54 (5.6%) [mism45 perr6] |
| ordering | 122/300 (40.7%) [mism30 unsup110 perr38] | 63/200 (31.5%) [mism79 unsup56 perr2] | 2/11 (18.2%) [mism6 unsup2 perr1] | 19/57 (33.3%) [mism18 unsup19 perr1] | 21/56 (37.5%) [mism17 unsup17 perr1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 389, 'exact': 1387, 'R2': 278, 'R3': 120}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 52 | 105 | 33 | 9 | 9 | 50 | 98 | 103 |
| e1 | 146 | 407 | 28 | 28 | 40 | 300 | 155 | 380 |
| e1-redteam | 0 | 21 | 3 | 1 | 2 | 8 | 3 | 5 |
| e2 | 28 | 43 | 13 | 1 | 17 | 18 | 47 | 29 |
| e4 | 14 | 58 | 17 | 1 | 14 | 13 | 47 | 36 |

biocpy heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 249/443 (56.2%) [mism73 unsup90 perr31] |
| empty | 125/223 (56.1%) [mism36 unsup42 perr20] |
| mcols | 447/804 (55.6%) [mism124 unsup187 perr46] |
| multi_seqlevel | 1006/1645 (61.2%) [mism229 unsup310 perr100] |
| na_arg | 128/222 (57.7%) [mism29 unsup53 perr12] |
| names | 366/694 (52.7%) [mism108 unsup166 perr54] |
| r_error | 217/440 (49.3%) [mism50 unsup173] |
| seqlengths | 338/593 (57.0%) [mism101 unsup121 perr33] |
| strand | 678/1187 (57.1%) [mism170 unsup265 perr74] |
| zero_width | 666/1210 (55.0%) [mism168 unsup306 perr70] |

## Runs

| port | tag | frozen sha (before/after) | clean | adapter sha256 | run s |
|---|---|---|---|---|---|
| e0 | frozen/e0-control | 18042e1/18042e1 | True/True | ca1bd8664566 | 1.9 |
| e0b | frozen/e0b-no-r | 3615722/3615722 | True/True | 5ec500b43049 | 4.1 |
| e1 | frozen/e1-differential | c5fa8df/c5fa8df | True/True | 94e31e16a6ae | 1.8 |
| e2 | frozen/e2-conformance | 70f0c0b/70f0c0b | True/True | d96f4c89f921 | 1.2 |
| e3 | frozen/e3-spec-adr | e0e4b41/e0e4b41 | True/True | 6b64d1b64459 | 1.4 |
| e4 | frozen/e4-primitive-cases | 9e2fa35/9e2fa35 | True/True | cb8663bc0e30 | 4.4 |
| biocpy | None | None/None | True/True | fe33e23a9756 | 40.6 |

## Self-tests

| adapter | scored | pass | mismatch | adapter_error | unsupported | port_error_unexpected |
|---|---|---|---|---|---|---|
| r_echo | 9650 | 9650 | 0 | 0 | 0 | 0 |
| mutant | 9650 | 7753 | 506 | 415 | 207 | 769 |

