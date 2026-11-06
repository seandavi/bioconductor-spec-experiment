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

## 2026-11-06 — E3 spec + ADRs done (harness v2)

| Stage | Wall time |
|---|---|
| plan | 2.9 min |
| implement (spec S1–S4 + code) | 45.6 min |
| review | 3.9 min |
| revise | 1.3 min |
| **total** | **54 min** |

- Spec: 85 clauses (later 86, GR.19 added in revise) in 11 files,
  ~4,400 words / 325 lines; 14 ADRs, ~1,400 words. ~2,020 lines src,
  ~1,780 lines tests, 177 tests, 21 commits. Verified: 177 passed;
  isolation 0 findings; reviewer wrote only its output.
- Traceability test: every live clause must be cited by a test and every
  cited id must exist (`@pytest.mark.spec(...)`). Allow-list empty.
- Handoff glitch: the planner wrote spec stages S1–S4 *for the
  implementer*, so the Sonnet implementer wrote the spec. In E3 the "spec
  author" is Sonnet, not Opus. (Not wrong per brief; note for comparison.)
- **Approach drift, third time: unprompted differential testing.** The
  implementer built a throwaway R fuzz (R generates random cases and
  answers, Python compared; ~110k comparisons) in `spec/scratch/`,
  *gitignored*. It found real corner cases now in spec + ADR-0013/14:
  `reduce(min_gapwidth=0)` drops repeated zero-width ranges; GRanges
  `intersect`/`setdiff` clip to seqlengths but `union` doesn't; per-op
  result-strand rules in parallel ops; `*` query nearest merges +/* and
  mirrored − candidates. The brief said a harness was "not required", not
  forbidden. Scripts (1,129 lines of R probes) copied to
  `ledger/evidence/e3-scratch/`; the 9.7 MB of TSVs left in the worktree.
- Consequence: E3's spec is better than "spec alone" would produce, but the
  evidence behind it is not in the repo. Its long-term artifact is prose +
  examples; the reproducible check was thrown away.
- Review: 1 major (masked/NA subscripts silently selected data, violating
  the spec's own CORE.9; an NA `find_overlaps(select="first")` result
  indexed row 0), 2 minor (spec vs tests disagreed on seqinfo merging in
  set ops: tests were defining undocumented semantics). All fixed; the
  spec gained GR.19.
- Pattern across E0, E3, E4: given R, agents *always* probe it heavily.
  The approaches differ less in whether R gets consulted than in **what
  evidence is kept**: nothing (E0 comments), prose examples (E3),
  re-runnable goldens (E4).

## 2026-11-06 — E2 formal conformance (with serialization) done (harness v2)

| Stage | Wall time |
|---|---|
| plan | 3.9 min |
| implement | 51.3 min |
| review | 3.2 min |
| revise | 3.3 min |
| **total** | **62 min** |

- ~2,270 lines package src; harness ~1,630 lines Python (codec, op
  adapters, runner, report) + ~430 lines R (codec, runner, selftest);
  FORMAT.md 574 words; 380 KB R-generated expected outputs; 921 cases over
  ~155 clauses; 34 commits. Verified: pytest 2037 passed, 4 xfailed;
  `Rscript conformance/r/run.R check` → 921 cases, 0 not passing.
  Isolation 0 findings; reviewer wrote only its output.
- **Serialization was cheap. My prediction was wrong.** Format: typed
  JSON (`$int/$dbl/$chr/$lgl` tags, NA as null, `$type` objects for
  IRanges/GRanges/Seqinfo/DFrame/Hits/Rle/RleList/IntegerList). R codec +
  JSON writer + selftest: ~1 commit, no redesign. Python codec: ~1 commit,
  797/797 goldens round-tripped on the first run. The agent's own
  verdict: "Biggest unexpected cost was not serialization but reading
  goldens to find which R behaviours are contracts and which are
  accidents" (~25 quirks).
- Where the real harness cost went: **per-op adapters** (`ops.py`, ~750
  lines, largest harness file). R calling conventions leak there: 1-based/
  negative/zero subscripts, NA arguments, `...` varargs, named vectors,
  `ignore.strand` only on GRanges, partial R defaults. Same layer E4 kept
  thin by `do.call` on the R side.
- Codec hacks that survived: `revmap` special-cased *by column name*
  (0-based list ↔ 1-based IntegerList); the reviewer found the general
  IntegerList/Factor mcols path broken (finding 3), fixed by typing columns
  independently of name.
- Committed its fuzz (`conformance/fuzz.py`, random R-vs-Python). Found 6
  behaviors the 900 hand-authored cases missed (reduce drop.empty.ranges,
  GRanges setdiff/intersect clipping, p* strand rules, nearest ties,
  precede/follow self). Fourth experiment to fuzz unprompted; the first
  to keep it. Overlap with E3's fuzz finds is high (clipping, p* strand),
  which suggests these are the real hotspots.
- Known gap: GRanges precede/follow/nearest with start/end ≤ 0 leak
  across sequences in R (sentinel number line); not cloned.
- Review: 3 major (mutable accessors: `x.width[0] = -5` worked;
  `c()` truncated double mcols to the first frame's int dtype; codec
  partial vs FORMAT.md), 2 minor (R runner ignored per-case compare mode;
  clause wording vs cases). All fixed, +16 cases.
- Mutable-accessor/aliasing bugs now seen in E0b, E2 (and E3 subscripts):
  the most common defect class across ports, and invisible to any
  R-vs-Python comparison. Only language-local review catches it.

## 2026-11-06 — E1 differential + red team done (harness v2)

| Stage | Model | Wall time |
|---|---|---|
| plan | Opus | 4.7 min |
| implement | Sonnet | 72.5 min |
| red team | Sol | 15.0 min |
| review | Sol | 5.4 min |
| revise | Sonnet | 35.6 min |
| **total** | | **133 min** |

- ~2,920 lines src; harness ~2,050 lines Python + ~205 lines R; 41k lines
  of case YAML (783 hand-written + 2,818 from 15 seeded generators, then
  red-team and revise additions); 4 MB committed R result cache (keyed by
  R/package versions + harness-script digest, so pytest runs without R);
  8 ADRs; 45 commits. Verified: 4,746 passed, 0 xfail (12 s, cache hits).
  Isolation: 1 audit hit, benign (a `cd` to a mistyped nonexistent dir,
  then into its own worktree). Red-team commits (00:42–00:52) touched only
  cases and tests/unit, never src/. Reviewer wrote only its output.
- **The red team was the highest-yield stage of the whole experiment.**
  After 3,601 differential cases were green, Sol wrote 119 adversarial
  cases in 15 min: 67 agreed, **52 diverged (47 bugs, 5 judgment calls)**.
  Clusters: circular seqinfo (17: overlaps across the origin, count/any/
  subset, nearest), metadata recycling/validation, empty/default inputs,
  `pintersect` existing `hit` column. Revise cleared all 52: fixes plus
  ADR-0007 (nearest family with coordinates ≤ 0) and ADR-0008 (circular
  maxgap).
- **Red team expanded scope.** E1 is the only port with circular-sequence
  semantics; E0, E0b and E3 list circular wrap-around as not done. TARGET
  lists `isCircular` but not wrap semantics explicitly; the red team
  treated R behavior as the scope. Scope = what the adversary tests.
- **Contract vs accident, live.** Reviewer suggested NA coverage weights
  stay local to their range. Revise rejected it with evidence: on 4,000
  random layouts, R has two internal engines; NA stays local only when the
  ranges tile the window, otherwise everything from the first NA range on
  is NA. E1 cloned that exactly (0 mismatches). This is Hervé's
  overspecification point in one case: R-normative-by-default copied what
  is almost certainly an implementation accident. A human would likely
  write an ADR saying "implementation-defined" instead.
- Cost: E1 is ~2× E2 and ~3× E0 in wall time. Most of the extra is
  implement (generators + R runs) and the revise round on 52 cases.
