# Protocol

Five independent experiments port the same target (TARGET.md) to Python
using different processes. We compare design, development velocity,
fidelity to R, and maintenance cost.

## Experiments

| Dir | Approach |
|---|---|
| `experiments/e0-control/` | Plain port with its own pytest suite. No extra process. |
| `experiments/e1-differential/` | Port with R normative; adversarial agents derive differential test cases (run in R and Python); divergences become discussion and ADRs. |
| `experiments/e2-conformance/` | Formal language-neutral conformance suite first, including serialization/deserialization of inputs and outputs; R and Python both run it. |
| `experiments/e3-spec-adr/` | Spec language and ADRs first. Implementation and tests are guided by, not defined by, the spec. What is not documented is not important. |
| `experiments/e4-primitive-cases/` | Shared case files with primitive inputs and projected primitive outputs (no object serialization); R runner fills goldens, Python runner checks them. |

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

## Evaluation (after all experiments finish)

1. Pool every divergence any experiment claims to have found.
2. Adjudicate each against actual R execution; record the verdict.
3. Score every Python implementation against the adjudicated pool
   (and BiocPy, as a bonus comparison).
4. Maintenance replay: apply a real historical upstream behavior change
   (from GenomicRanges/IRanges NEWS) and see how each approach detects
   and propagates it.
5. Qualitative: which artifacts would a Bioconductor developer read?
