# ADR-0005: coordinate overflow is an error, not NA coordinates

Status: accepted (divergence from R)

## Context
`shift(<GRanges>, 2147483647)` in R wraps to NA starts with only a warning and
returns a GRanges whose `start()` is NA (IRanges errors in the same situation).
Case: `gr-shift-overflow-error`.

## Options
1. Mirror R: return an object with NA coordinates and a warning.
2. Raise `ValueError` as IRanges does.

## Decision
Option 2. An object with NA coordinates violates the class invariant
(`start`/`width` are int64 without NA) and every downstream op would need NA
handling. IRanges and GRanges now behave the same.
