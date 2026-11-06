# Protocol

Five independent experiments port the same target (TARGET.md) to Python
using different processes. We compare design, development velocity,
fidelity to R, and maintenance cost.

## Experiments

| Dir | Approach |
|---|---|
| `experiments/e0-control/` | Plain port with its own pytest suite. No extra process. (Pilot; may consult R ad hoc.) |
| `experiments/e0b-no-r/` | Same as E0 but with no access to R: port from knowledge alone. |
| `experiments/e1-differential/` | Port with R normative; adversarial agents derive differential test cases (run in R and Python); divergences become discussion and ADRs. |
| `experiments/e2-conformance/` | Formal language-neutral conformance suite first, including serialization/deserialization of inputs and outputs; R and Python both run it. |
| `experiments/e3-spec-adr/` | Spec language and ADRs first. Implementation and tests are guided by, not defined by, the spec. What is not documented is not important. |
| `experiments/e4-primitive-cases/` | Shared case files with primitive inputs and projected primitive outputs (no object serialization); R runner fills goldens, Python runner checks them. |

## Pipeline

Every experiment runs `harness/run-experiment.js` with the same models:
plan (Opus 5.5) → implement (Sonnet 5.5) → [red team, E1 only (Sol)] →
review (Sol, no-edit, shell allowed) → revise (Sonnet, resumed). Only the
approach brief in `harness/approaches/` differs. E0 ran harness v1
(read-only reviewer without shell, no terseness rule); all later
experiments run v2.

## Isolation

- Each experiment lives on its own branch `exp/eN-*`, created from the
  `baseline` tag, in its own git worktree outside this directory.
- An experiment's agents read and write only inside their own
  experiment directory, plus TARGET.md and PROTOCOL.md. They do not read
  other branches, other worktrees, or `background/`.
- Each experiment is its own uv package (`pyproject.toml`, `src/`
  layout) with a distinct import name (`ranges_e0`, `ranges_e1`, ...).
- Same agent model for every experiment.

## History and velocity

- Commit at every meaningful stage: scaffold, each operation group,
  each test batch, each redesign. Never squash. Merges into `main` use
  `--no-ff`.
- Commit messages say what changed and *why*, especially for rewrites
  and redesigns ("rewrote overlap engine: interval tree → sorted sweep
  because ...").
- Velocity is measured, not budgeted: wall time, tokens, commits,
  rewrites, and any human intervention (logged in the experiment's
  `NOTES.md` with a timestamp).
- Safety ceiling: an experiment stops at its done criteria or after
  ~4 hours of agent wall time, whichever comes first. Stopping at the
  ceiling is a result, not a failure.

## Observations

The orchestrator keeps a running ledger on branch `ledger` (`LEDGER.md`):
observations, protocol changes and their reasons, and workflow ideas.

## Evaluation (after all experiments finish)

1. Pool every divergence any experiment claims to have found.
2. Adjudicate each against actual R execution; record the verdict.
3. Score every Python implementation against the adjudicated pool
   (and BiocPy, as a bonus comparison).
4. Maintenance replay: apply a real historical upstream behavior change
   (from GenomicRanges/IRanges NEWS) and see how each approach detects
   and propagates it.
5. Qualitative: which artifacts would a Bioconductor developer read?
