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

## Self-tests

| adapter | scored | pass | mismatch | adapter_error | unsupported | port_error_unexpected |
|---|---|---|---|---|---|---|
| r_echo | 9650 | 9650 | 0 | 0 | 0 | 0 |
| mutant | 9650 | 7753 | 506 | 415 | 207 | 769 |

