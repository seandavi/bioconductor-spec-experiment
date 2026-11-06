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
