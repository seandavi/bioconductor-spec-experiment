# ADR-0002: Python API surface

Status: accepted (plan stage)

## Context
R uses S4 generics with dotted argument names and copy-on-modify
replacement functions (`names(x) <- ...`).

## Options
1. Module-level functions mirroring R (`reduce(x)`).
2. Methods named after R generics in snake_case (`x.reduce(min_gapwidth=2)`).
3. Pythonic renames (`x.merge()`, `x.complement()`).

## Decision
Option 2. One-to-one name mapping keeps the differential adapter a
mechanical table and lets R users find things. Objects are immutable;
replacement functions become `set_*` methods returning copies. mcols are
`polars.DataFrame | None` (`None` = zero columns).
