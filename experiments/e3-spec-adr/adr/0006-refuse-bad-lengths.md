# ADR-0006 Refuse length mismatches instead of warn-and-recycle

Status: accepted · Clauses: CONV.6, CORE.3, SET.2, NEAR.1

## Context
R recycles mismatched vectors with a warning
(`IRanges(start=1:3, width=1:2)`), but errors for `pintersect` of length 2
vs 1. Warnings are easy to miss in Python.

## Options
1. Mimic R case by case.
2. Refuse every mismatch except length 1 (and `distance`, where R recycles by design).

## Decision
Option 2: only length-1 values recycle. Parallel set ops require equal length.

## Consequences
A few R calls that warn and succeed raise here. Spec states the rule once (CONV.6).
