# ADR-0003 mcols are polars DataFrames

Status: accepted · Clauses: CORE.8, CORE.12, INTER.5, INTER.8

## Context
R `DataFrame` carries row count even with zero columns and holds list
columns (`revmap`). A polars `DataFrame` with no columns has height 0.

## Options
1. pandas (heavier; not our preferred stack).
2. polars `DataFrame`; internal `None` for "no columns".
3. dict of numpy arrays (no list columns, no dtypes for NA).

## Decision
Option 2. `x.mcols` returns a `DataFrame` with `len(x)` rows when it has
columns, and an empty (zero-column) frame otherwise. `_mcols.py` is the only
place that handles the zero-column case. `revmap` is a `List(Int64)` column.

## Consequences
A zero-column mcols never carries a row count; `len(x)` is authoritative.
Concat uses a diagonal (union-of-columns) stack with nulls (CORE.12).
