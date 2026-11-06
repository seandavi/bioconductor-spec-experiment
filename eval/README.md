# Eval harness: fidelity to R

Scores Python ports of IRanges/GRanges against live R. Plan: `EVAL-PLAN.md`. Contract: `spec/`. Adapters: `ADAPTERS.md`. Deviations: `NOTES.md`.

## Run

| Step | Command | Output |
|---|---|---|
| convert E1/E2/E4 | `just convert` | `cases/{e1,e2,e4}.jsonl`, `unconvertible.jsonl`, `conversion-summary.json` |
| heldout | `just heldout` (then `just seeds`) | `cases/heldout.jsonl`, `heldout/COVERAGE.md` |
| R oracle | `just oracle` (~25 s cold, 4 workers) | `build/resolved.jsonl` (port-visible), `build/expected.jsonl` (comparator only) |
| self-tests | `just selftest` | r_echo 100 %, mutant == prediction, pytest |
| one port | `just run e1` | `build/results/e1.jsonl`, `results/scored/e1.jsonl` |
| report | `just report` | `results/report.md`, `results/matrix.json` |

## Layout

| Path | What |
|---|---|
| `src/bioc_eval/` | eval-venv tools: convert, oracle, compare, run, report, validate, seeds |
| `runner/` | `run_port.py`, `adapter_base.py`, `binding.py`: stdlib + numpy, run inside each port's venv |
| `oracle/oracle.R` | R side: eval decode, call, projection (own JSON writer, `%.17g`) |
| `adapters/` | `r_echo.py`, `mutant.py` now; one file per port next phase |
| `heldout/` | seeded generator (stdlib), `COVERAGE.md`, `seeds.txt` |
| `ports/biocpy/` | pinned uv env: iranges 0.7.3, genomicranges 0.9.0 |
| `selftest/fakeport/` | stdlib package that only raises |
| `tests/` | unit + end-to-end pipeline tests (`just test`) |

Port results are never scored without `ports.toml`'s frozen sha and a clean worktree before and after.
