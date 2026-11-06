# Ledger

Orchestrator's running record: observations, protocol changes, evidence.
Newest entries at the bottom. Times are local (America/Denver).
Workflow and agent-loop ideas live in [WORKFLOW-IDEAS.md](WORKFLOW-IDEAS.md).

## 2026-11-05 — setup

- Repo initialised. `baseline` tag = TARGET.md + PROTOCOL.md only.
  `background/` (Zulip threads, ChatGPT chat) committed after the tag, so
  experiment branches never contain it (keeps the seqlevel-order hint out).
- Pipeline fixed for all experiments: plan Opus 5.5 → implement Sonnet 5.5
  → review Sol (gpt-6.1-sol via `openai-codex`) → revise (Sonnet, resumed).
- Codex CLI with a ChatGPT login rejects `gpt-6.1-sol`; the pi
  `openai-codex` provider (after Sean authenticated) serves it. Built-in
  `codex-exec` agent can't take a model override, so it would have run the
  Codex default (gpt-5.6-terra, medium).

## 2026-11-05 — E0 pilot (harness v1, d72e8bd)

| Stage | Model | Wall time |
|---|---|---|
| plan | Opus 5.5 | 3.7 min |
| implement | Sonnet 5.5 | 27.5 min |
| review | Sol | 3.6 min (incl. two supervisor round-trips) |
| revise | Sonnet | 1.3 min |
| **total** | | **36 min** |

- Output: ~2,650 lines src, ~2,270 lines tests, 20 commits, 306 tests
  passing (re-run by me). Every TARGET.md group implemented; some
  non-required extras (distanceToNearest, circular wrap) skipped.
- Isolation held: parsed all tool calls in the four transcripts; no path
  outside the E0 worktree, no `background/`, no BiocPy, no other refs.
- **The control did informal differential testing on its own.** ~280 R
  invocations, a reusable `r.sh` helper, 110 test assertions annotated
  `# R: ...`. Allowed by the brief ("ad hoc R snippets"), but it means E0
  measures "agent + R REPL", not "agent alone". Added E0b (no R) to cover
  the knowledge-only case.
- Review value: Sol found 4 issues with Python-only probes against E0's
  *own* stated contracts (coverage weight recycling crashes, distance
  silently recycles length 2 vs 3, flank vector `start` gives width of
  wrong length, nearest self-match on empty range). 300 green tests missed
  all four.
- Revise rejected finding 1 with R evidence (R also returns the empty range
  as its own nearest, and duplicates the hit). I confirmed in R:
  `nearest(IRanges(5,width=0))` → `1`; `select="all"` → `(1,1)` twice.
  A miniature ADR loop: dispute settled by running the reference.
- Protocol issue: the built-in `reviewer` agent has no shell, so Sol asked
  me to run pytest and probes (two supervisor requests). Fixed in v2.

## 2026-11-05 — harness v2 (70ee971); E0b, E1–E4 launched in parallel

- Changes: reviewer = no-edit `worker` with shell (audit `git status`
  after); red-team stage for E1 (Sol writes differential cases meant to
  break the port, failing ones → fix or ADR); `rAccess=false` for E0b;
  terse-writing rule for all agent prose (Sean's request).
- Confounds to remember when comparing to E0: v2 terseness rule (affects
  prose volume), reviewer can now probe on its own (may find more), E1 has
  an extra adversarial stage (more agent time by design).
- E1–E4 briefs are short on purpose; planners choose harness/format/spec
  design (Sean: "allow the planner to do its thing").
- Run ids: e0b 756556f2, e1 447d7fce, e2 206d6ce4, e3 7b277843, e4 01fcc0b8.

## 2026-11-05 — E0b no-R done (harness v2)

| Stage | Wall time |
|---|---|
| plan | 4.2 min |
| implement | 16.1 min |
| review | 3.5 min |
| revise | 1.4 min |
| **total** | **25 min** |

- ~2,070 lines src, ~1,850 lines tests, 24 commits, 211 tests passing
  (re-run by me). Worktree clean.
- Isolation audit with `--no-r`: 0 findings. No R use at all. Reviewer's
  only writes were its own output artifact.
- Faster than E0 (16 vs 27.5 min implement): no R round trips. Smaller
  suite (211 vs 306), but the planner wrote down the R rules from memory
  for every operation (PLAN §4) and tagged doubtful ones `UNVERIFIED`;
  21 such tags in tests/src/NOTES. Listed risks: zero-width overlaps,
  `resize(fix="center")` rounding, `narrow(width=)`, GRanges `gaps` with NA
  seqlengths, strand in `type=start/end` and precede/follow, `pintersect`
  on adjacent ranges, `within` + `maxgap`.
- Review: 3 major + 1 minor, all fixed, none rejected. The majors were
  design and completeness issues (no `ignore_strand` on resize/flank/
  promoters, which TARGET requires; mutable shared mcols despite a
  documented immutable design; `with_seqinfo` couldn't drop unused levels),
  not R semantics.
- Early read: without R the agent becomes *explicit about its uncertainty*.
  The UNVERIFIED list is a ready-made target list for differential testing.
  Whether E0b is actually less faithful than E0 waits for the pooled
  evaluation.

## 2026-11-06 — E4 primitive cases done (harness v2)

| Stage | Wall time |
|---|---|
| plan | 4.9 min |
| implement | 33.5 min |
| review | 3.6 min |
| revise | 3.6 min |
| **total** | **46 min** |

- ~1,900 lines src, ~1,230 lines tests, ~220 lines of runners (R + Python),
  ~4,000 lines of case YAML, 240 KB of R-generated JSON goldens.
  715 cases, 305 clauses (`clauses.md`), 11 ADRs, 29 commits.
- Verified by me: `uv run pytest` 2308 passed, 6 skipped;
  `Rscript runners/r/run_cases.R --check` → goldens up to date (this is
  the upstream-drift detector: rerun after an R upgrade). Isolation audit
  0 findings; reviewer wrote only /tmp probe + its output.
- The case format stayed primitive, but **"no serialization" still had a
  cross-language parsing gotcha**: R's yaml reads unquoted `y`/`n`/`Y`
  as logicals, PyYAML doesn't (an input named `y` became TRUE). Also ints
  beyond int32 become NA in R yaml. Fixed with yaml handlers + a quoting
  rule. Small, but even primitives need a written convention.
- R runner is ~generic: `do.call(op, args)` with R names verbatim, no op
  table in R. The op mapping lives only on the Python side.
- ADRs (11) are mostly API conventions (0-based indices, masked-array NA,
  no negative subscripts, errors compared by class, recycle-from-one) plus
  real semantic calls: `oneof` comparator for `select="arbitrary"` and
  nearest ties (ADR-0005), and **R defects not ported** (ADR-0010: gaps/
  setdiff/disjoin at end = 2^31−1, maxgap = 2^31−1, zero-length
  precede/follow returning Hits). E4 found R edge-case defects.
- Independent confirmation of the Zulip seqlevel claim: E4's planner,
  probing R, recorded "GRanges seqlevels follow order of appearance".
  That matches what Sean's ADR example said about R (not yet checked
  against BiocPy here).
- Review: 4 major + 3 minor, all fixed. Best finding was about the
  *harness*, not the port: "exact" comparison applied float tolerance to
  integers, so an off-by-one coordinate would pass. A conformance suite's
  comparator is itself code that needs review and tests.
