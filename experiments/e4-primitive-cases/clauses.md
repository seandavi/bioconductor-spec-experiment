# Clauses

One behavior per row. Append-only; retire with `(retired)`.
Every clause needs >=1 case (checked by `tests/test_case_inventory.py`).

| id | statement | R ref |
|---|---|---|
| IR-CON-START-WIDTH-1 | `IRanges(start, width)`: `end = start + width - 1` | `?IRanges` |
| IR-CON-START-END-1 | `IRanges(start, end)`: `width = end - start + 1` | `?IRanges` |
| IR-CON-END-WIDTH-1 | `IRanges(end, width)`: `start = end - width + 1` | `?IRanges` |
| IR-CON-EMPTY-1 | zero-length IRanges is valid; all accessors return empty vectors | `?IRanges` |
| IR-CON-ZEROWIDTH-1 | zero-width ranges are valid: `end = start - 1` | `?IRanges` |
| IR-CON-NEGWIDTH-1 | negative width is an error | `?IRanges` |
| IR-CON-NAMES-1 | names are carried through; NA names allowed | `?IRanges` |
| IR-CON-MCOLS-1 | mcols columns are carried through, in order | `?mcols` |
