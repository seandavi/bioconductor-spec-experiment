# Heldout generator brief

Write `eval/heldout/generate.py` (stdlib only) emitting
`eval/cases/heldout.jsonl` in the eval case schema. Scope comes from
`TARGET.md` **only**.

## Isolation (hard)

- Read: `TARGET.md`, this file, `spec/CASES.md`, `spec/case.schema.json`, `spec/ops.json`.
- Do not read: `spec/CONVERTERS.md`, `eval/cases/` (other than your output), any experiment worktree or branch, any port's tests/cases/goldens, `background/`, R output of any kind.
- The generator never runs R and never looks at results. R decides every expected value later.

## Seeds

```python
seed = int.from_bytes(hashlib.sha256(f"bioc-spec-eval/heldout/v1/{op}/{k}".encode()).digest()[:8], "big")
rng = random.Random(seed)
```

64-bit, derived from a namespace no experiment uses. The build stage greps
every experiment worktree for each seed value and fails on any hit.

## Volume and shape

- Every op in `ops.json`, each in IRanges **and** GRanges form where the op applies to both (`Seqinfo` ops: Seqinfo and GRanges).
- ~30 cases per (op, class); ~60 for `findOverlaps`, `nearest`, `reduce`, `coverage`. Target 2-3k cases.
- Ids `heldout/<op>/<class>/<k>`; `source: heldout`; `origin: {suite: heldout, file: "heldout/generate.py", id: "<op>/<k>", tags: [features...]}`.
- Inputs small (n = 0..12) on a narrow coordinate window (≈ -5..60) so overlaps, adjacency and nesting occur; ~2 % of cases near int32 limits.
- Valid inputs (R must build them): non-negative widths, strand in `+ - *`, seqnames within seqinfo.
- Kwargs: sample R's legal arguments (`ops.json` `named`) with realistic values; ~10 % deliberately invalid values (negative widths, unknown enum strings, wrong lengths) so error agreement is measured.

## Required coverage (assert in the generator; print a table)

Per op family, at least 5 cases each with:

| Feature | Meaning |
|---|---|
| `empty` | a range input of length 0 |
| `zero_width` | a width-0 range |
| `names` | named ranges (incl. duplicated names) |
| `mcols` | 1-2 metadata columns (int, dbl, chr, lgl; some NA) |
| `strand` | GRanges with mixed `+ - *`; and `ignore.strand` where the op takes it |
| `seqinfo` | multiple seqlevels, unused levels, non-sorted level order, some known seqlengths, a circular level, a genome |
| `na_arg` | an NA in an atomic argument |

Plus: constructor ops (any two of start/end/width; names; error combos),
setters (`names<-`, `mcols<-`, `seqlevels<-`, `seqlengths<-`, ...), `c` of
inputs with differing seqinfo (merge), `trim` against seqlengths, ordering
across seqlevels and strand.

## Done

Deterministic (two runs byte-identical); every line validates against
`case.schema.json`; coverage table printed and committed as
`eval/heldout/COVERAGE.md`.
