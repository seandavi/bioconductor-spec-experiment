# bioconductor-spec-experiment

How much process does a faithful Python port of a Bioconductor class
need? Six independent agentic ports of IRanges/GRanges, each built with a
different process, scored against live R on 9,650 cases, then put through
a simulated upstream change.

- [TARGET.md](TARGET.md): what gets built, and when it counts as done
- [PROTOCOL.md](PROTOCOL.md): the experiments, isolation rules, metrics
- [QUESTIONS.md](QUESTIONS.md): R behaviors that look like accidents, where a person should decide what a port must match
- `background/`: a paraphrased summary of the Bioconductor Zulip discussion that prompted this

## Experiments

Same pipeline for all: Opus 5.5 plans → Sonnet 5.5 implements → Sol
(gpt-6.1-sol) reviews → Sonnet revises. Only the process brief differs
(`harness/approaches/`).

| | Process | Build time | Adjusted held-out fidelity | Detects upstream change |
|---|---|---|---|---|
| E0 | plain port; ad hoc R allowed | 36 min | 96.0% | no |
| E0b | plain port; no R | 25 min | 92.4% | no |
| E1 | differential harness + Sol red team | 133 min | **98.7%** | yes (version-keyed) |
| E2 | formal conformance suite first (typed JSON) | 62 min | 96.9% | yes |
| E3 | spec + ADRs first | 54 min | 94.2% | no |
| E4 | shared primitive cases, R goldens | 46 min | 96.2% | yes |
| BiocPy 0.9.0 | existing community port | — | 81.3% | — |

Fidelity = pass / (supported − R accidents − declared divergences) on 3,510
held-out cases no port had seen. Raw pass rates, per-source and per-family
tables: `eval/results/report.md`; adjudication: `eval/adjudication/ADJUDICATION.md`.

## Findings in one screen

- Fidelity tracks how much *recorded, adversarial* checking against R a
  process did. E1's red team found 52 divergences after 3,601 green cases.
- Given R, agents always probe it. Processes differ in what evidence they
  keep: nothing (E0), prose (E3), re-runnable goldens (E2, E4, E1).
- Serialization (E2) was cheap; translating R calling conventions and
  deciding contract-vs-accident were the real costs.
- Maintenance: given a release note, every approach propagated a
  `distance()` change in minutes with no regressions. Only E1/E2/E4
  *noticed* it without being told. Detection is where process pays.
- BiocPy sorts seqlevels lexically; R and all six ports use first
  appearance (the Zulip thread's example, confirmed).
- "R is normative" copies R's accidents along with its contracts. Six
  places where that needs a human call are in [QUESTIONS.md](QUESTIONS.md).

## Layout

| Path | What |
|---|---|
| `experiments/<exp>/` | each port as evaluated (tags `frozen/<exp>`) |
| `eval/` | fidelity harness, adapters, oracle, results, adjudication; `eval/replay/` replay scoring |
| `harness/` | workflow scripts, approach briefs, isolation audit, replay patch |
| `ledger/` | orchestrator's running record, review copies, workflow ideas |

To re-run the evaluation, recreate the worktrees it expects with
`sh harness/recreate-worktrees.sh`, then `cd eval && just oracle && just selftest && just run e1`
(`just oracle` regenerates the gitignored `build/` from live R in about 25 s).
R 4.6 with Bioc 3.24 and uv are required; paths in `eval/ports.toml` are absolute.

Branches `exp/*`, `replay/*`, `eval*`, `ledger` are merged with `--no-ff`.
Replay merges were reverted on `main` (their ports implement the simulated
change, not real R); the updated ports live on `replay/<exp>`.
