# ranges_e4

Independent Python port of Bioconductor `IRanges` / `GRanges` (scope: `TARGET.md`).
Experiment **e4**: conformance is checked by **shared case files of
primitives**. R is normative. Objects are never serialized.

```
uv run pytest        # Python suite: unit tests + one test per conformance case
just goldens         # R regenerates goldens/ from cases/
just check-goldens   # R recomputes and diffs against committed goldens (exit 1 on drift)
```

Needs `Rscript` with IRanges, GenomicRanges, yaml, jsonlite (R 4.6 / Bioc 3.24
when written) for the two `just` recipes only. `uv run pytest` needs no R.

## How conformance works

```
cases/<iranges|granges>/<group>.yaml      hand-authored: primitive inputs + call + compare modes
        |                  \
        | Rscript            \  Python runner (tests/conformance/)
        v                     v
goldens/.../<group>.json     build objects natively -> run op -> project -> compare
(R builds objects natively,   ^
 runs, projects; never edited)|__ same projection as R (project.R <-> project.py)
```

- **Case** (`id`, `clause`, `inputs`, `call`, optional `compare`, `oneof_from`,
  `adr`, `python`): inputs are constructor arguments of ints, strings,
  logicals (`IRanges{start,end,width,names,mcols}`,
  `GRanges{seqnames,start,...,strand,seqinfo}`, `Seqinfo{...}`, `int|num|chr|lgl`
  literals). `call` = R op name + positional input names + R-named kwargs.
  **No expected values in case files.**
- **Golden**: R's result projected onto primitive vectors by class
  (`IRanges`: start/end/width/names/mcols; `GRanges`: + seqnames/strand/seqinfo;
  `Hits`; `Rle`; `RleList`; `Seqinfo`; atomic `vector`; or `error`).
  R warnings are recorded, not compared.
- **Compare modes** (default `exact`): `ignore`, `oneof` (R computes allowed
  candidates by re-running with `select="all"`), `multiset` (rows), `set`.
  Any non-exact mode, and any `python:` override (`expect_error`, `skip`),
  must cite an ADR file in `adr/` (enforced by `tests/test_case_inventory.py`).
- **Clauses** (`clauses.md`): one behavior per row, id `IR-/GR-/SI-<GROUP>-<OP>-<n>`.
  Each case names a clause; each clause needs a case.
- **Index base**: Python indices are 0-based, goldens stay R-native. The Python
  runner adds 1 when projecting index fields (`ops.py` `index` flag, `Hits`,
  `revmap`) and subtracts 1 from `int` subscript inputs (ADR-0001).

Size at close-out: 701 cases (389 IRanges, 312 GRanges), 305 clauses, 11 ADRs,
2250 passing tests, 6 skipped (ADR-0010).

## Package design

Seam: every algorithm is a kernel over plain `start`/`width` int64 arrays;
`GRanges` only decides *which* subset of ranges a kernel sees and how to
reassemble the result.

| Module | Role |
|---|---|
| `iranges.py` | `IRanges` (start, width, names, polars `mcols`): construct, subset, concat, public ops |
| `granges.py` | `GRanges` (seqlevel codes, strand codes, `IRanges`, mcols, `Seqinfo`) |
| `seqinfo.py` | `Seqinfo` + `merge` (NA fills, conflicts raise) |
| `intra.py` | shift/narrow/resize/flank/promoters kernels; a `reverse` mask gives strand awareness |
| `inter.py` | range/reduce/gaps/disjoin: one sorted sweep (`start - run_end - 1 < min_gapwidth`) |
| `setops.py` | union/intersect/setdiff via coverage breakpoints; parallel `p*` ops |
| `overlaps.py`, `nearest.py` | predicates pinned to R on grids; sorted-prefix candidate bound |
| `_group.py` | run a kernel per (seqlevel, strand) group (`ignore_strand` -> `*`) |
| `coverage.py`, `rle.py` | difference-array coverage (circular wrap), minimal `Rle` |
| `ordering.py` | stable lexsort keys: IRanges (start,width); GRanges (seqlevel,strand,start,width) |

API: methods on the classes, snake_case (`find_overlaps`, `subset_by_overlaps`,
`is_circular`); R `.` -> `_` in argument names (`min_gapwidth`, `ignore_strand`).
`IRanges.concat(*xs)` / `GRanges.concat(*xs)`. Named R vectors (`seqlengths`,
`isCircular`, `genome`) are `dict`s. NA integer results are `np.ma.MaskedArray`.
Deliberate divergences: ADR-0001 (0-based), 0004, 0006, 0007, 0008, 0009, 0010, 0011.

## Maintaining the artifacts

| Change | Steps |
|---|---|
| New behavior or op | add clause to `clauses.md` -> add case(s) -> `just goldens` -> implement until `uv run pytest` is green |
| R upgrade / upstream behavior change | `just check-goldens` lists drifted files. Re-run `just goldens`, commit goldens **alone**, read the diff like code, then fix Python or add an ADR |
| Python intentionally differs | keep R's golden; add `python: {expect_error: true}` or `{skip: ...}` + `adr:` + ADR file. Never edit a golden |
| Implementation-defined result | `compare` mode + ADR (see ADR-0005) |

Authoring rules:
- Cases hold inputs only; goldens are generated, never hand-edited.
- Quote strings that YAML 1.1 reads as booleans (`yes`, `no`, `on`, `off`).
  `y`/`n` are kept as strings by the R runner's handlers (PyYAML never read them as booleans).
- Keep coordinates within int32 (R's yaml turns bigger ints into NA).
- Avoid index `0` and negative subscripts in cases (ADR-0004).
- Probe R before writing a clause; note surprising R output in `NOTES.md`.

`tests/test_case_inventory.py` guards: unique ids, clause <-> case <-> golden
coverage, ADR citations, every TARGET operation has a clause and cases in both
classes. Beyond the cases, `tests/test_*.py` hold brute-force oracles on seeded
random inputs (reduce, disjoin, gaps, set algebra, overlaps, nearest, coverage,
ordering) and Python-only behavior (`test_edge.py`).

## Known limits

- `GRanges` `nearest`/`precede`/`follow` need an explicit subject (R's no-subject form is not ported).
- Candidate search is O(n*m) worst case per overlap/nearest call (correctness first).
- `Rle` is only the coverage result; `RleList` is a `dict[str, Rle]`.
- Not ported: R exclusion subscripts (ADR-0004), cyclic recycling (ADR-0007), lone-`start` `IRanges` form (ADR-0011).
