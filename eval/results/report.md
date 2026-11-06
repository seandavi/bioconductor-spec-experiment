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
| e0 | 2773/3510 (79.0%) [mism149 unsup567 perr21] | 3180/4385 (72.5%) [mism579 unsup500 perr126] | 69/119 (58.0%) [mism27 unsup9 perr14] | 836/921 (90.8%) [mism34 unsup36 perr15] | 672/715 (94.0%) [mism30 unsup8 perr5] | 7530/9650 (78.0%) |
| e0b | 2505/3510 (71.4%) [mism174 unsup755 perr76] | 2606/4385 (59.4%) [mism627 unsup691 perr461] | 60/119 (50.4%) [mism23 unsup21 perr15] | 765/921 (83.1%) [mism58 unsup54 perr44] | 603/715 (84.3%) [mism68 unsup13 perr31] | 6539/9650 (67.8%) |
| e1 | 2984/3510 (85.0%) [mism75 unsup444 perr7] | *3940/4385 (89.9%) [mism370 unsup73 perr2] | *110/119 (92.4%) [mism9] | 877/921 (95.2%) [mism18 unsup26] | 706/715 (98.7%) [mism6 perr3] | 8617/9650 (89.3%) |
| e2 | 2891/3510 (82.4%) [mism116 unsup487 perr16] | 3128/4385 (71.3%) [mism649 unsup239 perr369] | 84/119 (70.6%) [mism27 unsup4 perr4] | *891/921 (96.7%) [mism15 unsup15] | 685/715 (95.8%) [mism23 unsup4 perr3] | 7679/9650 (79.6%) |
| e3 | 2571/3510 (73.2%) [mism142 unsup745 perr52] | 2916/4385 (66.5%) [mism532 unsup700 perr237] | 61/119 (51.3%) [mism25 unsup18 perr15] | 805/921 (87.4%) [mism31 unsup64 perr21] | 639/715 (89.4%) [mism39 unsup21 perr16] | 6992/9650 (72.5%) |
| e4 | 2642/3510 (75.3%) [mism121 unsup717 perr30] | 2950/4385 (67.3%) [mism885 unsup382 perr168] | 61/119 (51.3%) [mism30 unsup14 perr14] | 815/921 (88.5%) [mism30 unsup62 perr14] | *698/715 (97.6%) [mism8 unsup1 perr8] | 7166/9650 (74.3%) |
| biocpy | 2192/3510 (62.5%) [mism398 unsup796 perr124] | 1944/4385 (44.3%) [mism1339 unsup761 perr341] | 41/119 (34.5%) [mism38 unsup18 perr22] | 570/921 (61.9%) [mism172 unsup152 perr27] | 419/715 (58.6%) [mism161 unsup86 perr49] | 5166/9650 (53.5%) |

Cell: `pass/scored (rate) [mism = mismatch, adapt = adapter_error, unsup = unsupported, perr = port_error_unexpected]`.

## e0: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 221/240 (92.1%) [mism9 unsup7 perr3] | 280/399 (70.2%) [mism19 unsup31 perr69] | 3/10 (30.0%) [mism4 unsup1 perr2] | 50/60 (83.3%) [mism3 unsup2 perr5] | 39/46 (84.8%) [mism5 unsup2] |
| accessors | 409/480 (85.2%) [mism4 unsup67] | 32/36 (88.9%) [mism3 unsup1] | - | 25/27 (92.6%) [mism2] | 14/15 (93.3%) [unsup1] |
| metadata | 104/120 (86.7%) [mism11 perr5] | 21/24 (87.5%) [mism2 unsup1] | 9/21 (42.9%) [perr12] | 9/11 (81.8%) [perr2] | 2/2 (100.0%) |
| seqinfo | 315/570 (55.3%) [unsup255] | 107/114 (93.9%) [unsup7] | 6/8 (75.0%) [mism2] | 30/50 (60.0%) [mism1 unsup19] | 25/25 (100.0%) |
| subset_concat | 77/120 (64.2%) [mism1 unsup42] | 351/462 (76.0%) [mism33 unsup78] | 1/1 (100.0%) | 50/60 (83.3%) [unsup10] | 36/37 (97.3%) [unsup1] |
| intra | 281/300 (93.7%) [mism4 unsup3 perr12] | 321/363 (88.4%) [mism14 perr28] | 2/2 (100.0%) | 124/126 (98.4%) [perr2] | 82/88 (93.2%) [mism5 perr1] |
| inter | 240/300 (80.0%) [mism14 unsup46] | 336/405 (83.0%) [mism58 unsup7 perr4] | 12/12 (100.0%) | 91/101 (90.1%) [mism6 perr4] | 101/104 (97.1%) [mism3] |
| setops | 338/360 (93.9%) [mism16 unsup6] | 397/441 (90.0%) [mism36 perr8] | 4/4 (100.0%) | 102/107 (95.3%) [mism5] | 106/108 (98.1%) [mism1 perr1] |
| overlaps | 273/300 (91.0%) [mism27] | 539/848 (63.6%) [mism304 unsup3 perr2] | 8/23 (34.8%) [mism15] | 185/195 (94.9%) [mism10] | 109/113 (96.5%) [mism3 unsup1] |
| nearest | 256/300 (85.3%) [mism43 perr1] | 528/582 (90.7%) [mism49 perr5] | 10/15 (66.7%) [mism5] | 72/76 (94.7%) [mism3 perr1] | 62/67 (92.5%) [mism4 perr1] |
| coverage | 68/120 (56.7%) [mism6 unsup46] | 105/511 (20.5%) [mism26 unsup372 perr8] | 3/12 (25.0%) [mism1 unsup8] | 43/51 (84.3%) [mism3 unsup4 perr1] | 41/54 (75.9%) [mism8 unsup3 perr2] |
| ordering | 191/300 (63.7%) [mism14 unsup95] | 163/200 (81.5%) [mism35 perr2] | 11/11 (100.0%) | 55/57 (96.5%) [mism1 unsup1] | 55/56 (98.2%) [mism1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 255, 'exact': 327, 'R2': 197, 'R3': 40}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 3 | 24 | 22 | 10 | 0 | 62 | 11 | 17 |
| e1 | 1 | 238 | 19 | 26 | 0 | 169 | 1 | 125 |
| e1-redteam | 0 | 17 | 4 | 1 | 0 | 6 | 0 | 2 |
| e2 | 2 | 2 | 10 | 0 | 1 | 8 | 3 | 8 |
| e4 | 2 | 10 | 2 | 0 | 0 | 10 | 1 | 5 |

e0 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 321/443 (72.5%) [mism32 unsup87 perr3] |
| empty | 185/223 (83.0%) [mism13 unsup24 perr1] |
| mcols | 626/804 (77.9%) [mism47 unsup130 perr1] |
| multi_seqlevel | 1202/1645 (73.1%) [mism102 unsup333 perr8] |
| na_arg | 50/222 (22.5%) [mism53 unsup111 perr8] |
| names | 538/694 (77.5%) [mism48 unsup107 perr1] |
| r_error | 282/440 (64.1%) [mism62 unsup96] |
| seqlengths | 408/593 (68.8%) [mism43 unsup141 perr1] |
| strand | 882/1187 (74.3%) [mism81 unsup217 perr7] |
| zero_width | 969/1210 (80.1%) [mism60 unsup176 perr5] |

## e0b: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 215/240 (89.6%) [mism9 unsup13 perr3] | 273/399 (68.4%) [mism15 unsup31 perr80] | 6/10 (60.0%) [mism1 unsup1 perr2] | 48/60 (80.0%) [mism2 unsup5 perr5] | 40/46 (87.0%) [mism3 unsup2 perr1] |
| accessors | 404/480 (84.2%) [mism9 unsup67] | 31/36 (86.1%) [mism5] | - | 24/27 (88.9%) [mism3] | 15/15 (100.0%) |
| metadata | 115/120 (95.8%) [perr5] | 24/24 (100.0%) | 9/21 (42.9%) [perr12] | 9/11 (81.8%) [perr2] | 2/2 (100.0%) |
| seqinfo | 315/570 (55.3%) [unsup255] | 107/114 (93.9%) [unsup7] | 6/8 (75.0%) [mism2] | 30/50 (60.0%) [mism1 unsup19] | 25/25 (100.0%) |
| subset_concat | 56/120 (46.7%) [unsup45 perr19] | 225/462 (48.7%) [mism1 unsup78 perr158] | 0/1 (0.0%) [mism1] | 37/60 (61.7%) [unsup11 perr12] | 33/37 (89.2%) [mism1 perr3] |
| intra | 228/300 (76.0%) [mism6 unsup63 perr3] | 301/363 (82.9%) [mism15 unsup19 perr28] | 0/2 (0.0%) [unsup2] | 122/126 (96.8%) [unsup2 perr2] | 79/88 (89.8%) [mism6 unsup2 perr1] |
| inter | 185/300 (61.7%) [mism33 unsup79 perr3] | 196/405 (48.4%) [mism122 unsup80 perr7] | 9/12 (75.0%) [unsup3] | 80/101 (79.2%) [mism13 unsup5 perr3] | 87/104 (83.7%) [mism12 unsup5] |
| setops | 282/360 (78.3%) [mism44 unsup13 perr21] | 321/441 (72.8%) [mism34 unsup10 perr76] | 2/4 (50.0%) [unsup2] | 78/107 (72.9%) [mism18 unsup1 perr10] | 82/108 (75.9%) [mism19 perr7] |
| overlaps | 275/300 (91.7%) [mism15 unsup10] | 501/848 (59.1%) [mism287 unsup58 perr2] | 8/23 (34.8%) [mism13 unsup2] | 187/195 (95.9%) [mism8] | 108/113 (95.6%) [mism4 unsup1] |
| nearest | 231/300 (77.0%) [mism52 perr17] | 396/582 (68.0%) [mism108 perr78] | 10/15 (66.7%) [mism5] | 60/76 (78.9%) [mism10 perr6] | 43/67 (64.2%) [mism11 perr13] |
| coverage | 71/120 (59.2%) [mism3 unsup46] | 99/511 (19.4%) [mism22 unsup366 perr24] | 3/12 (25.0%) [mism1 unsup8] | 43/51 (84.3%) [mism2 unsup3 perr3] | 38/54 (70.4%) [mism11 unsup1 perr4] |
| ordering | 128/300 (42.7%) [mism3 unsup164 perr5] | 132/200 (66.0%) [mism18 unsup42 perr8] | 7/11 (63.6%) [unsup3 perr1] | 47/57 (82.5%) [mism1 unsup8 perr1] | 51/56 (91.1%) [mism1 unsup2 perr2] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 182, 'exact': 501, 'R2': 221, 'R3': 46}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 7 | 65 | 18 | 12 | 7 | 42 | 0 | 27 |
| e1 | 2 | 266 | 15 | 30 | 58 | 117 | 1 | 159 |
| e1-redteam | 0 | 15 | 3 | 1 | 0 | 5 | 0 | 2 |
| e2 | 3 | 22 | 9 | 0 | 5 | 6 | 1 | 13 |
| e4 | 2 | 36 | 1 | 0 | 5 | 12 | 1 | 15 |

e0b heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 289/443 (65.2%) [mism31 unsup111 perr12] |
| empty | 170/223 (76.2%) [mism12 unsup35 perr6] |
| mcols | 542/804 (67.4%) [mism57 unsup173 perr32] |
| multi_seqlevel | 1076/1645 (65.4%) [mism100 unsup429 perr40] |
| na_arg | 80/222 (36.0%) [mism28 unsup111 perr3] |
| names | 480/694 (69.2%) [mism54 unsup135 perr25] |
| r_error | 281/440 (63.9%) [mism42 unsup117] |
| seqlengths | 377/593 (63.6%) [mism32 unsup171 perr13] |
| strand | 786/1187 (66.2%) [mism83 unsup287 perr31] |
| zero_width | 847/1210 (70.0%) [mism81 unsup243 perr39] |

## e1: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 223/240 (92.9%) [mism10 unsup6 perr1] | 399/399 (100.0%) | 10/10 (100.0%) | 57/60 (95.0%) [mism3] | 46/46 (100.0%) |
| accessors | 409/480 (85.2%) [mism4 unsup67] | 32/36 (88.9%) [mism4] | - | 25/27 (92.6%) [mism2] | 15/15 (100.0%) |
| metadata | 120/120 (100.0%) | 24/24 (100.0%) | 21/21 (100.0%) | 11/11 (100.0%) | 2/2 (100.0%) |
| seqinfo | 416/570 (73.0%) [mism1 unsup149 perr4] | 114/114 (100.0%) | 8/8 (100.0%) | 35/50 (70.0%) [unsup15] | 25/25 (100.0%) |
| subset_concat | 77/120 (64.2%) [mism1 unsup42] | 388/462 (84.0%) [mism1 unsup73] | 1/1 (100.0%) | 50/60 (83.3%) [unsup10] | 37/37 (100.0%) |
| intra | 298/300 (99.3%) [mism1 perr1] | 362/363 (99.7%) [perr1] | 2/2 (100.0%) | 126/126 (100.0%) | 85/88 (96.6%) [perr3] |
| inter | 261/300 (87.0%) [unsup39] | 405/405 (100.0%) | 12/12 (100.0%) | 101/101 (100.0%) | 102/104 (98.1%) [mism2] |
| setops | 336/360 (93.3%) [mism18 unsup6] | 441/441 (100.0%) | 4/4 (100.0%) | 102/107 (95.3%) [mism5] | 107/108 (99.1%) [mism1] |
| overlaps | 294/300 (98.0%) [mism6] | 835/848 (98.5%) [mism13] | 22/23 (95.7%) [mism1] | 188/195 (96.4%) [mism7] | 111/113 (98.2%) [mism2] |
| nearest | 272/300 (90.7%) [mism28] | 582/582 (100.0%) | 12/15 (80.0%) [mism3] | 76/76 (100.0%) | 67/67 (100.0%) |
| coverage | 76/120 (63.3%) [mism3 unsup40 perr1] | 159/511 (31.1%) [mism351 perr1] | 7/12 (58.3%) [mism5] | 51/51 (100.0%) | 54/54 (100.0%) |
| ordering | 202/300 (67.3%) [mism3 unsup95] | 199/200 (99.5%) [mism1] | 11/11 (100.0%) | 55/57 (96.5%) [mism1 unsup1] | 55/56 (98.2%) [mism1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 11, 'R2': 14, 'exact': 443, 'R3': 10}

| source | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|
| heldout | 23 | 19 | 10 | 1 | 3 | 3 | 18 |
| e1 | 3 | 13 | 0 | 0 | 3 | 180 | 341 |
| e1-redteam | 2 | 1 | 0 | 0 | 1 | 3 | 5 |
| e2 | 1 | 9 | 0 | 2 | 2 | 0 | 4 |
| e4 | 2 | 1 | 0 | 0 | 2 | 0 | 1 |

e1 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 368/443 (83.1%) [mism17 unsup57 perr1] |
| empty | 201/223 (90.1%) [mism3 unsup18 perr1] |
| mcols | 670/804 (83.3%) [mism24 unsup109 perr1] |
| multi_seqlevel | 1359/1645 (82.6%) [mism48 unsup232 perr6] |
| na_arg | 157/222 (70.7%) [mism3 unsup60 perr2] |
| names | 576/694 (83.0%) [mism37 unsup80 perr1] |
| r_error | 358/440 (81.4%) [mism3 unsup79] |
| seqlengths | 472/593 (79.6%) [mism19 unsup101 perr1] |
| strand | 1001/1187 (84.3%) [mism42 unsup141 perr3] |
| zero_width | 1042/1210 (86.1%) [mism30 unsup136 perr2] |

## e2: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 231/240 (96.2%) [unsup6 perr3] | 337/399 (84.5%) [mism27 perr35] | 6/10 (60.0%) [mism3 perr1] | 59/60 (98.3%) [unsup1] | 44/46 (95.7%) [mism1 perr1] |
| accessors | 404/480 (84.2%) [mism9 unsup67] | 31/36 (86.1%) [mism5] | - | 24/27 (88.9%) [mism3] | 15/15 (100.0%) |
| metadata | 120/120 (100.0%) | 24/24 (100.0%) | 21/21 (100.0%) | 11/11 (100.0%) | 2/2 (100.0%) |
| seqinfo | 449/570 (78.8%) [mism5 unsup113 perr3] | 113/114 (99.1%) [mism1] | 7/8 (87.5%) [mism1] | 46/50 (92.0%) [unsup4] | 25/25 (100.0%) |
| subset_concat | 74/120 (61.7%) [mism1 unsup45] | 351/462 (76.0%) [mism33 unsup78] | 0/1 (0.0%) [mism1] | 51/60 (85.0%) [unsup9] | 36/37 (97.3%) [mism1] |
| intra | 239/300 (79.7%) [mism10 unsup44 perr7] | 308/363 (84.8%) [mism33 unsup14 perr8] | 2/2 (100.0%) | 126/126 (100.0%) | 85/88 (96.6%) [mism2 perr1] |
| inter | 242/300 (80.7%) [mism11 unsup47] | 303/405 (74.8%) [mism60 unsup42] | 12/12 (100.0%) | 101/101 (100.0%) | 102/104 (98.1%) [mism2] |
| setops | 334/360 (92.8%) [mism13 unsup13] | 388/441 (88.0%) [mism44 unsup1 perr8] | 3/4 (75.0%) [mism1] | 107/107 (100.0%) | 106/108 (98.1%) [mism1 perr1] |
| overlaps | 271/300 (90.3%) [mism13 unsup16] | 494/848 (58.3%) [mism252 unsup100 perr2] | 8/23 (34.8%) [mism13 unsup2] | 188/195 (96.4%) [mism7] | 107/113 (94.7%) [mism2 unsup4] |
| nearest | 256/300 (85.3%) [mism44] | 508/582 (87.3%) [mism49 perr25] | 10/15 (66.7%) [mism5] | 73/76 (96.1%) [mism3] | 63/67 (94.0%) [mism4] |
| coverage | 74/120 (61.7%) [mism2 unsup41 perr3] | 124/511 (24.3%) [mism94 unsup4 perr289] | 4/12 (33.3%) [mism3 unsup2 perr3] | 50/51 (98.0%) [mism1] | 45/54 (83.3%) [mism9] |
| ordering | 197/300 (65.7%) [mism8 unsup95] | 147/200 (73.5%) [mism51 perr2] | 11/11 (100.0%) | 55/57 (96.5%) [mism1 unsup1] | 55/56 (98.2%) [mism1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 293, 'exact': 337, 'R2': 156, 'R3': 44}

| source | kind | length | names | not_in_universe | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|
| heldout | 8 | 7 | 18 | 11 | 57 | 1 | 14 |
| e1 | 2 | 191 | 15 | 29 | 222 | 1 | 190 |
| e1-redteam | 0 | 17 | 3 | 1 | 6 | 1 | 3 |
| e2 | 3 | 0 | 9 | 0 | 2 | 1 | 0 |
| e4 | 2 | 9 | 1 | 0 | 6 | 1 | 4 |

e2 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 361/443 (81.5%) [mism22 unsup58 perr2] |
| empty | 190/223 (85.2%) [mism9 unsup24] |
| mcols | 643/804 (80.0%) [mism26 unsup134 perr1] |
| multi_seqlevel | 1325/1645 (80.5%) [mism73 unsup239 perr8] |
| na_arg | 119/222 (53.6%) [mism32 unsup61 perr10] |
| names | 553/694 (79.7%) [mism43 unsup97 perr1] |
| r_error | 305/440 (69.3%) [mism57 unsup78] |
| seqlengths | 455/593 (76.7%) [mism30 unsup107 perr1] |
| strand | 933/1187 (78.6%) [mism59 unsup188 perr7] |
| zero_width | 992/1210 (82.0%) [mism51 unsup166 perr1] |

## e3: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 210/240 (87.5%) [mism9 unsup12 perr9] | 290/399 (72.7%) [mism24 perr85] | 4/10 (40.0%) [mism3 unsup1 perr2] | 47/60 (78.3%) [mism2 unsup5 perr6] | 38/46 (82.6%) [mism3 perr5] |
| accessors | 401/480 (83.5%) [mism12 unsup67] | 28/36 (77.8%) [mism8] | - | 22/27 (81.5%) [mism4 unsup1] | 14/15 (93.3%) [mism1] |
| metadata | 115/120 (95.8%) [perr5] | 24/24 (100.0%) | 9/21 (42.9%) [perr12] | 9/11 (81.8%) [perr2] | 2/2 (100.0%) |
| seqinfo | 314/570 (55.1%) [unsup256] | 107/114 (93.9%) [unsup7] | 5/8 (62.5%) [unsup3] | 29/50 (58.0%) [unsup21] | 25/25 (100.0%) |
| subset_concat | 74/120 (61.7%) [mism1 unsup45] | 351/462 (76.0%) [mism33 unsup78] | 0/1 (0.0%) [mism1] | 49/60 (81.7%) [unsup11] | 36/37 (97.3%) [mism1] |
| intra | 213/300 (71.0%) [mism6 unsup72 perr9] | 272/363 (74.9%) [mism35 unsup21 perr35] | 0/2 (0.0%) [mism2] | 120/126 (95.2%) [unsup4 perr2] | 78/88 (88.6%) [mism5 unsup4 perr1] |
| inter | 232/300 (77.3%) [mism11 unsup50 perr7] | 313/405 (77.3%) [mism38 unsup43 perr11] | 12/12 (100.0%) | 96/101 (95.0%) [unsup1 perr4] | 101/104 (97.1%) [mism2 unsup1] |
| setops | 303/360 (84.2%) [mism41 unsup13 perr3] | 368/441 (83.4%) [mism54 unsup10 perr9] | 2/4 (50.0%) [unsup2] | 94/107 (87.9%) [mism11 unsup1 perr1] | 95/108 (88.0%) [mism12 perr1] |
| overlaps | 273/300 (91.0%) [mism11 unsup16] | 496/848 (58.5%) [mism250 unsup100 perr2] | 8/23 (34.8%) [mism13 unsup2] | 188/195 (96.4%) [mism7] | 107/113 (94.7%) [mism2 unsup4] |
| nearest | 239/300 (79.7%) [mism47 perr14] | 457/582 (78.5%) [mism49 perr76] | 10/15 (66.7%) [mism5] | 67/76 (88.2%) [mism6 perr3] | 56/67 (83.6%) [mism5 perr6] |
| coverage | 43/120 (35.8%) [unsup77] | 68/511 (13.3%) [mism14 unsup418 perr11] | 3/12 (25.0%) [mism1 unsup8] | 35/51 (68.6%) [unsup14 perr2] | 34/54 (63.0%) [mism7 unsup12 perr1] |
| ordering | 154/300 (51.3%) [mism4 unsup137 perr5] | 142/200 (71.0%) [mism27 unsup23 perr8] | 8/11 (72.7%) [unsup2 perr1] | 49/57 (86.0%) [mism1 unsup6 perr1] | 53/56 (94.6%) [mism1 perr2] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 227, 'exact': 336, 'R2': 160, 'R3': 46}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 7 | 41 | 21 | 12 | 0 | 46 | 1 | 14 |
| e1 | 2 | 197 | 18 | 30 | 0 | 161 | 1 | 123 |
| e1-redteam | 0 | 15 | 3 | 1 | 0 | 7 | 0 | 2 |
| e2 | 3 | 11 | 10 | 0 | 1 | 3 | 1 | 2 |
| e4 | 2 | 19 | 2 | 0 | 0 | 10 | 1 | 5 |

e3 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 296/443 (66.8%) [mism28 unsup110 perr9] |
| empty | 175/223 (78.5%) [mism11 unsup35 perr2] |
| mcols | 590/804 (73.4%) [mism39 unsup166 perr9] |
| multi_seqlevel | 1109/1645 (67.4%) [mism93 unsup423 perr20] |
| na_arg | 70/222 (31.5%) [mism29 unsup109 perr14] |
| names | 493/694 (71.0%) [mism55 unsup138 perr8] |
| r_error | 273/440 (62.0%) [mism46 unsup121] |
| seqlengths | 387/593 (65.3%) [mism32 unsup170 perr4] |
| strand | 809/1187 (68.2%) [mism74 unsup286 perr18] |
| zero_width | 884/1210 (73.1%) [mism69 unsup242 perr15] |

## e4: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 213/240 (88.8%) [mism9 unsup12 perr6] | 299/399 (74.9%) [mism19 perr81] | 4/10 (40.0%) [mism3 unsup1 perr2] | 49/60 (81.7%) [mism1 unsup5 perr5] | 45/46 (97.8%) [perr1] |
| accessors | 404/480 (84.2%) [mism9 unsup67] | 31/36 (86.1%) [mism5] | - | 24/27 (88.9%) [mism3] | 15/15 (100.0%) |
| metadata | 85/120 (70.8%) [unsup30 perr5] | 21/24 (87.5%) [unsup3] | 9/21 (42.9%) [perr12] | 8/11 (72.7%) [unsup1 perr2] | 2/2 (100.0%) |
| seqinfo | 314/570 (55.1%) [unsup256] | 107/114 (93.9%) [unsup7] | 5/8 (62.5%) [unsup3] | 29/50 (58.0%) [unsup21] | 25/25 (100.0%) |
| subset_concat | 74/120 (61.7%) [mism1 unsup45] | 351/462 (76.0%) [mism33 unsup78] | 1/1 (100.0%) | 49/60 (81.7%) [unsup11] | 37/37 (100.0%) |
| intra | 263/300 (87.7%) [mism4 unsup22 perr11] | 287/363 (79.1%) [mism27 unsup13 perr36] | 0/2 (0.0%) [mism2] | 124/126 (98.4%) [perr2] | 84/88 (95.5%) [perr4] |
| inter | 233/300 (77.7%) [mism19 unsup48] | 256/405 (63.2%) [mism75 unsup72 perr2] | 10/12 (83.3%) [mism2] | 90/101 (89.1%) [mism7 unsup4] | 102/104 (98.1%) [mism2] |
| setops | 317/360 (88.1%) [mism23 unsup13 perr7] | 369/441 (83.7%) [mism46 unsup18 perr8] | 2/4 (50.0%) [unsup2] | 97/107 (90.7%) [mism6 unsup1 perr3] | 107/108 (99.1%) [mism1] |
| overlaps | 277/300 (92.3%) [mism11 unsup12] | 501/848 (59.1%) [mism268 unsup77 perr2] | 8/23 (34.8%) [mism13 unsup2] | 188/195 (96.4%) [mism7] | 111/113 (98.2%) [mism2] |
| nearest | 235/300 (78.3%) [mism39 unsup25 perr1] | 447/582 (76.8%) [mism39 unsup70 perr26] | 8/15 (53.3%) [mism5 unsup2] | 58/76 (76.3%) [mism5 unsup12 perr1] | 63/67 (94.0%) [mism3 perr1] |
| coverage | 72/120 (60.0%) [mism5 unsup43] | 149/511 (29.2%) [mism347 unsup4 perr11] | 5/12 (41.7%) [mism5 unsup2] | 49/51 (96.1%) [mism1 perr1] | 52/54 (96.3%) [perr2] |
| ordering | 155/300 (51.7%) [mism1 unsup144] | 132/200 (66.0%) [mism26 unsup40 perr2] | 9/11 (81.8%) [unsup2] | 50/57 (87.7%) [unsup7] | 55/56 (98.2%) [unsup1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 192, 'exact': 670, 'R2': 174, 'R3': 38}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 8 | 35 | 18 | 8 | 0 | 37 | 1 | 14 |
| e1 | 2 | 453 | 15 | 26 | 8 | 147 | 1 | 244 |
| e1-redteam | 0 | 21 | 3 | 1 | 0 | 6 | 1 | 2 |
| e2 | 3 | 6 | 9 | 0 | 1 | 1 | 1 | 9 |
| e4 | 2 | 2 | 1 | 0 | 0 | 1 | 1 | 1 |

e4 heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 300/443 (67.7%) [mism25 unsup114 perr4] |
| empty | 176/223 (78.9%) [mism11 unsup33 perr3] |
| mcols | 604/804 (75.1%) [mism34 unsup163 perr3] |
| multi_seqlevel | 1129/1645 (68.6%) [mism75 unsup433 perr8] |
| na_arg | 69/222 (31.1%) [mism36 unsup110 perr7] |
| names | 504/694 (72.6%) [mism50 unsup137 perr3] |
| r_error | 282/440 (64.1%) [mism37 unsup121] |
| seqlengths | 389/593 (65.6%) [mism27 unsup175 perr2] |
| strand | 817/1187 (68.8%) [mism63 unsup301 perr6] |
| zero_width | 918/1210 (75.9%) [mism54 unsup230 perr8] |

## biocpy: family x source

| family | heldout | e1 | e1-redteam | e2 | e4 |
|---|---|---|---|---|---|
| construct | 145/240 (60.4%) [mism17 unsup61 perr17] | 242/399 (60.7%) [mism24 unsup40 perr93] | 6/10 (60.0%) [mism2 perr2] | 36/60 (60.0%) [mism3 unsup13 perr8] | 27/46 (58.7%) [mism6 unsup7 perr6] |
| accessors | 401/480 (83.5%) [mism12 unsup67] | 28/36 (77.8%) [mism8] | - | 23/27 (85.2%) [mism4] | 13/15 (86.7%) [mism2] |
| metadata | 104/120 (86.7%) [mism6 perr10] | 20/24 (83.3%) [mism4] | 4/21 (19.0%) [perr17] | 9/11 (81.8%) [perr2] | 1/2 (50.0%) [mism1] |
| seqinfo | 467/570 (81.9%) [mism29 unsup62 perr12] | 37/114 (32.5%) [mism76 unsup1] | 5/8 (62.5%) [mism3] | 30/50 (60.0%) [mism13 unsup6 perr1] | 24/25 (96.0%) [mism1] |
| subset_concat | 52/120 (43.3%) [mism16 unsup45 perr7] | 303/462 (65.6%) [mism68 unsup78 perr13] | 0/1 (0.0%) [mism1] | 41/60 (68.3%) [mism8 unsup7 perr4] | 24/37 (64.9%) [mism11 perr2] |
| intra | 216/300 (72.0%) [mism7 unsup76 perr1] | 303/363 (83.5%) [mism35 unsup25] | 0/2 (0.0%) [mism2] | 117/126 (92.9%) [mism5 unsup4] | 73/88 (83.0%) [mism7 unsup2 perr6] |
| inter | 242/300 (80.7%) [mism17 unsup27 perr14] | 302/405 (74.6%) [mism90 perr13] | 8/12 (66.7%) [mism3 perr1] | 93/101 (92.1%) [mism6 perr2] | 93/104 (89.4%) [mism5 perr6] |
| setops | 113/360 (31.4%) [mism62 unsup157 perr28] | 108/441 (24.5%) [mism154 unsup163 perr16] | 1/4 (25.0%) [mism1 unsup2] | 45/107 (42.1%) [mism16 unsup41 perr5] | 36/108 (33.3%) [mism22 unsup36 perr14] |
| overlaps | 119/300 (39.7%) [mism65 unsup107 perr9] | 299/848 (35.3%) [mism383 unsup158 perr8] | 7/23 (30.4%) [mism10 unsup6] | 111/195 (56.9%) [mism44 unsup40] | 70/113 (61.9%) [mism27 unsup15 perr1] |
| nearest | 187/300 (62.3%) [mism49 unsup42 perr22] | 217/582 (37.3%) [mism215 unsup65 perr85] | 8/15 (53.3%) [mism5 unsup2] | 43/76 (56.6%) [mism16 unsup13 perr4] | 36/67 (53.7%) [mism17 unsup1 perr13] |
| coverage | 6/120 (5.0%) [mism69 unsup42 perr3] | 22/511 (4.3%) [mism203 unsup175 perr111] | 0/12 (0.0%) [mism5 unsup6 perr1] | 3/51 (5.9%) [mism39 unsup9] | 1/54 (1.9%) [mism45 unsup8] |
| ordering | 140/300 (46.7%) [mism49 unsup110 perr1] | 63/200 (31.5%) [mism79 unsup56 perr2] | 2/11 (18.2%) [mism6 unsup2 perr1] | 19/57 (33.3%) [mism18 unsup19 perr1] | 21/56 (37.5%) [mism17 unsup17 perr1] |

Mismatch kinds (a case counts once per kind): rule used {'row4': 388, 'exact': 1363, 'R2': 278, 'R3': 79}

| source | kind | length | names | not_in_universe | order | r_error_port_value | type | value |
|---|---|---|---|---|---|---|---|---|
| heldout | 13 | 111 | 32 | 9 | 9 | 49 | 96 | 115 |
| e1 | 26 | 406 | 28 | 52 | 40 | 300 | 154 | 452 |
| e1-redteam | 0 | 21 | 3 | 1 | 2 | 8 | 3 | 5 |
| e2 | 4 | 43 | 13 | 2 | 17 | 18 | 47 | 43 |
| e4 | 3 | 58 | 17 | 2 | 14 | 13 | 47 | 39 |

biocpy heldout by input feature:

| feature | pass/scored |
|---|---|
| circular | 252/443 (56.9%) [mism70 unsup97 perr24] |
| empty | 125/223 (56.1%) [mism35 unsup46 perr17] |
| mcols | 452/804 (56.2%) [mism121 unsup199 perr32] |
| multi_seqlevel | 1028/1645 (62.5%) [mism214 unsup332 perr71] |
| na_arg | 127/222 (57.2%) [mism29 unsup56 perr10] |
| names | 379/694 (54.6%) [mism100 unsup177 perr38] |
| r_error | 206/440 (46.8%) [mism49 unsup185] |
| seqlengths | 342/593 (57.7%) [mism98 unsup128 perr25] |
| strand | 704/1187 (59.3%) [mism152 unsup285 perr46] |
| zero_width | 681/1210 (56.3%) [mism151 unsup326 perr52] |

## Runs

| port | tag | frozen sha (before/after) | clean | adapter sha256 | run s |
|---|---|---|---|---|---|
| e0 | frozen/e0-control | 18042e1/18042e1 | True/True | 8a8dec79a236 | 3.6 |
| e0b | frozen/e0b-no-r | 3615722/3615722 | True/True | 7c767f2273df | 5.6 |
| e1 | frozen/e1-differential | c5fa8df/c5fa8df | True/True | 28c01ec4ff7c | 2.8 |
| e2 | frozen/e2-conformance | 70f0c0b/70f0c0b | True/True | 7e5b545dba05 | 2.5 |
| e3 | frozen/e3-spec-adr | e0e4b41/e0e4b41 | True/True | 5f4b93fb9311 | 2.7 |
| e4 | frozen/e4-primitive-cases | 9e2fa35/9e2fa35 | True/True | 8ef3dd27f66b | 5.7 |
| biocpy | None | None/None | True/True | 70e3a7630a45 | 40.5 |

## Self-tests

| adapter | scored | pass | mismatch | adapter_error | unsupported | port_error_unexpected |
|---|---|---|---|---|---|---|
| r_echo | 9650 | 9650 | 0 | 0 | 0 | 0 |
| mutant | 9650 | 7753 | 506 | 415 | 207 | 769 |

