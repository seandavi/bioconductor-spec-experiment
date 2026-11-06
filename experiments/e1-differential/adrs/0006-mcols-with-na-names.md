# ADR-0006: `mcols()` works when names contain NA

Status: accepted (divergence from R)

## Context
R's `mcols(x)` builds a DataFrame whose rownames are `names(x)`, so it fails
("missing values not allowed in rownames") if any name is NA. Everything else
(print, subsetting, `c()`, ops) works on such objects. Case: `ir-mcols-accessor-na-names`.

## Options
1. Raise in `mcols` for NA names, as R.
2. Return the columns; Python mcols are polars frames without rownames.

## Decision
Option 2. Python mcols carry no rownames, so there is nothing to fail on. The
differential encoder reads R's mcols with `use.names=FALSE` so objects with NA names
still compare; only the explicit accessor call is an accepted divergence.
