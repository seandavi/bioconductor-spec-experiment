# Eval plan: scoring six ports (+ BiocPy) against live R

Goal: score E0, E0b, E1, E2, E3, E4 and BiocPy on fidelity to R (IRanges
2.47.x / GenomicRanges 1.65.x, R 4.6). R is normative; expected values come
only from the eval's own R runner.

## Contract (written in this stage)

| Artifact | Defines |
|---|---|
| `spec/CASES.md` + `case.schema.json` | case = inputs + call, R vocabulary; literal decode rules |
| `spec/PROJECTION.md` + `projection.schema.json` | one normalized shape per result kind (R conventions) |
| `spec/COMPARATOR.md` | outcomes, exact rule, the only three relaxations |
| `spec/ops.json` | op vocabulary, families, R formals (from R `formals()`), index-valued parts |
| `ADAPTERS.md` | adapter interface, binding, allowed vs forbidden, classification, self-tests |
| `spec/CONVERTERS.md` | E1/E2/E4 → eval mapping, unconvertible reasons |
| `spec/HELDOUT.md` | heldout generator brief (TARGET-only) |
| `ports.toml` | port registry: venv project, import, frozen tag/sha |

## Architecture

```
E1 / E2 / E4 suites ──convert──▶ cases/{e1,e2,e4}.jsonl     cases/heldout.jsonl ◀── generate (TARGET only)
                                  └──────────────┬──────────────┘
                                   just oracle   │  Rscript, one batch, cache per (versions, script, case hash)
                       ┌─────────────────────────┴───────────────────────┐
              build/resolved.jsonl                              build/expected.jsonl
              canonical inputs (R-built, projected) + calls      R result projection, universe, r_class, warnings
                       │  (the only file a port sees)                    │
   just run PORT: uv run --frozen --project <port> python runner/run_port.py --port PORT
                       │                                                 │
              build/results/PORT.jsonl ──────── compare (eval venv) ─────┘
                                                      │
                                    results/scored/PORT.jsonl ── just report ──▶ results/report.md, matrix.json
```

Key choices:

| Choice | Why |
|---|---|
| Ports get **R-canonical inputs** (R builds every input, projects it; adapters build from the projection) | an input-decoding quirk can't pollute every case; construction is tested where it belongs |
| Construction tested via constructor **ops** (`IRanges(...)`, `GRanges(...)`); converters rewrite source `identity(literal)` cases into constructor trees | keeps E1/E4 construct + bad-arg cases |
| Port runner stdlib + numpy + port only; JSONL between stages | runs unchanged in 7 different venvs |
| Expected values never reach the port venv | `resolved.jsonl` has none; only `r_echo`/`mutant` read `expected.jsonl` |
| Base binds R positional/named args via `ops.json` before the adapter sees them | one arg-matching implementation, not seven |
| Outcome classification by traceback (did the exception enter `port_packages`?) | mechanical split of `adapter_error` vs port error |
| Doubles exact, oracle writes `%.17g` | `jsonlite digits=NA` is 15 digits (verified); ulp-only diffs flagged `float_ulp`, not forgiven |
| Headline = heldout | E1/E2/E4 ports have home advantage on their own suites |

## Layout

```
eval/
  EVAL-PLAN.md ADAPTERS.md ports.toml justfile pyproject.toml   # eval tools venv (pyyaml, jsonschema, pytest)
  spec/                         contract (this stage)
  src/bioc_eval/                convert/{e1,e2,e4}.py, oracle.py, compare.py, report.py, run.py
  oracle/oracle.R               jsonlite + IRanges + GenomicRanges only
  runner/run_port.py, adapter_base.py     stdlib + numpy; executed inside each port venv
  adapters/r_echo.py mutant.py  (+ 7 port adapters, next phase)
  selftest/fakeport/            stdlib package that only raises
  heldout/generate.py COVERAGE.md
  ports/biocpy/                 uv project: iranges==0.7.3, genomicranges==0.9.0
  cases/                        committed pool, unconvertible.jsonl, conversion-summary.json
  build/                        gitignored: oracle cache, resolved/expected, raw results
  results/                      committed at the end: oracle snapshot, scored/*.jsonl, report.md, matrix.json
```

## Stages (each ends in a commit; never squash/amend/rebase)

| # | Stage | Deliverable | Accept when |
|---|---|---|---|
| S0 | Design | this file + contract artifacts + schema examples | examples validate against both schemas |
| S1 | Scaffold | `pyproject.toml` (uv, src layout), `justfile` (all recipes; unbuilt ones exit 1 with a message), `.gitignore`, `just validate` | `just validate` green on `spec/examples` |
| S2 | Converters | `convert/{e1,e2,e4}.py` + tests; `cases/` | every source case is converted or listed; read = converted + unconvertible per source; all lines schema-valid |
| S3 | Oracle | `oracle/oracle.R`, `oracle.py` (batching, cache, resolved/expected split, `domain`/`features`) | each case `ok` or `excluded` with reason; two runs byte-identical; projections schema-valid; R3 `universe` and R2 `sorted_by_query` present where required; 10 hand-checked cases |
| S4 | Heldout | **fresh agent** with `spec/HELDOUT.md` only (the designer has read the suites) | deterministic; coverage asserts pass; seed grep over all experiment worktrees finds nothing |
| S5 | Runner + comparator | `run_port.py`, `adapter_base.py`, `r_echo`, `mutant`, `fakeport`, `compare.py`, unit tests (binding, classification, each rule) | `just selftest`: r_echo 100 % pass; mutant outcome == predicted on every case |
| S6 | Report | `report.py`: matrix port x family x source x outcome; per-source headline; mismatch kinds; features; exclusions | counts reconcile (scored + excluded + unconvertible = read) on r_echo/mutant |
| S7 | BiocPy project | `ports/biocpy` with lockfile | `uv run --frozen` imports both pins |
| S8 | Adapters (next phase) | 7 fresh agents, one port each, from `ADAPTERS.md` + public API | per-adapter done list in ADAPTERS.md; review gate audits forbidden patterns |
| S9 | Full run | `just oracle && just run PORT` x 7 `&& just report` | report committed with versions, frozen SHAs, adapter digests |

Oracle notes (S3): one `Rscript` per batch; `options(warn = 1)`; warnings
recorded with `withCallingHandlers` and **not** muffled; each case in a fresh
environment; `universe` = same call with `select = "all"`; cache key =
sha256(R + IRanges + GenomicRanges + S4Vectors + Seqinfo versions,
`oracle.R` digest, case hash).

Run notes (S5/S9): `uv run --frozen` (never touches `uv.lock`); before and
after each run assert port worktree HEAD == `ports.toml` sha and clean
`git status`; record both in `results/PORT.meta.json`. Results append per
case with flush; `--resume` skips finished ids; per-case timeout 30 s.

## Justfile (S1)

| Recipe | Does |
|---|---|
| `just convert` | E1/E2/E4 → `cases/` |
| `just heldout` | regenerate `cases/heldout.jsonl` |
| `just validate` | schema-check examples, cases, resolved, expected |
| `just oracle` | R on all cases (cached) → `build/resolved.jsonl`, `build/expected.jsonl` |
| `just run PORT` | verify frozen sha → run port in its venv → score |
| `just roundtrip PORT` / `just adapter-check PORT` | adapter dev loop (no R values shown) |
| `just selftest` | r_echo + mutant + unit tests |
| `just report` | matrix + per-case JSONL + report.md |

## Risks

| Risk | Mitigation |
|---|---|
| Adapter quality confounds port score | thin-adapter rules; round-trip check; traceback classification; `adapter_error` reported separately; review gate |
| Adapter "fixes" port output | forbidden table in ADAPTERS.md; adapters never see R values in the dev loop; reviewer diff audit |
| Adapter wrongly feeds the port → port raises → counted as port error | round trip + `adapter-check` surfaces `port_error_unexpected` concentrations per op for adapter review |
| Converter changes a source case's meaning | converters mirror each source's own R decoder explicitly (CONVERTERS.md); unconvertible listed with reasons |
| E1 dominates the pool (4.5k of ~6k converted) | never a pooled-only score: per source columns, heldout headline |
| Home-suite advantage (E1 port on e1 cases, ...) | report marks the diagonal; heldout is blind for all |
| Heldout contaminated by suite knowledge | fresh agent, isolation list, seed namespace + grep |
| Hidden R nondeterminism | double oracle run must be identical |
| Hits order relaxation masks a real bug | R2 keeps query order exact and is keyed on R's class, not op name; mutant control case |
| Float ulp noise counted as mismatch | kept exact by contract; `float_ulp` kind lets the report separate it |
| Slow ports (O(n·m)) hang on large cases | per-case timeout → `port_error_unexpected` (`origin: timeout`) |
| `uv run` mutating a frozen port | `--frozen`; sha + clean check before/after |
| Input-kind gaps (port lacks Factor/IntegerList mcols) | `Unsupported("input:<kind>")`, visible per feature |

## Decisions to confirm (non-blocking; current choice in bold)

- Timeout: **folded into `port_error_unexpected`** (category list is fixed) vs a sixth category.
- `unsupported` in the pass-rate denominator: **yes** (no incentive to dodge).
- Port error while building an R-valid input: **always `port_error_unexpected`**, never an R1 pass.
- Factor inside Rle (`seqnames()`, `strand()`) projects as **chr vector**; plain factors keep levels.
