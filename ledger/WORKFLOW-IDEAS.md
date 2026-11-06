# Workflow and agent-loop ideas

Ideas that came up while running the experiments. Each says what was seen
and what to try. Not yet adopted unless marked.

## 1. "Read-only" should mean no edits, not no shell — adopted in v2

Seen: the built-in read-only reviewer couldn't run tests, so it asked the
supervisor to run pytest and four probe scripts. Each round trip stalls the
child and costs orchestrator context.
Try: reviewers run as normal workers told not to edit; the orchestrator
checks `git status` afterwards. Gives the reviewer real evidence and keeps
the "no writes" guarantee checkable.

## 2. Controls are defined by capability, not by instruction

Seen: E0 was allowed ad hoc R and used it ~280 times, quietly turning the
control into informal differential testing. An agent uses every tool it
has.
Try: define experimental arms by what tools exist (PATH without Rscript, a
sandbox profile), not by what the brief discourages. Where that's not
possible (E0b), audit transcripts after the fact.

## 3. Audit isolation from transcripts, automatically

Seen: a 30-line transcript parse (tool-call paths outside the worktree,
forbidden packages, `git log --all`) answered "did isolation hold?" in
seconds.
Try: make it a harness step that runs after every child and fails the run
on violations. Cheaper than sandboxing and leaves evidence.

## 4. Reject review findings only with evidence from the reference

Seen: E0's revise rejected a finding by running R and recording the output.
That's the right shape: a disagreement settled by the reference rather than
by argument.
Try: revise prompts require reproducible evidence (R command and output)
for every rejected finding. Rejections without it count as unresolved.

## 5. Probe against the code's own contracts

Seen: Sol's best findings came from probing E0's *own* docstrings and
README claims (recycling rules, lengths), not from R. 300 tests missed them
because the implementer's tests share its blind spots.
Try: a cheap "contract probe" stage: a different model reads public
docstrings only and writes probes for the stated contracts. Language-local
and needs no R.

## 6. Hand reviewers the facts up front

Seen: the reviewer's first supervisor request was for test results and
commit history the harness already had.
Try: the review brief includes the implementer's final reply (tests,
commits). The reviewer re-runs rather than discovers.

## 7. Persist artifacts out of retention-managed storage

Seen: review outputs and receipts live in temp/retention-managed
directories.
Try: the orchestrator copies them into the repo (here, `ledger/reviews/`,
`ledger/runs/`) as soon as a stage finishes.

## 8. Ask implementers to tag their uncertainty, then aim tests at it

Seen: E0b (no R) tagged 21 expectations `UNVERIFIED` and listed its riskiest
semantics. E0 (with R) checked ad hoc and left no such list.
Try: every implement stage must emit a machine-readable uncertainty list
(operation, argument combo, what it assumed). A differential or red-team
stage takes that list as its first targets. Cheap, and it points adversarial
effort where the author's knowledge is weakest.

## 9. Review the oracle, not just the subject

Seen: E4's worst defect was in its comparator (integer "exact" match with
float tolerance). Every green case was weaker than it looked.
Try: a review checklist item for any experiment with a harness: mutate
one golden / one Python result and confirm the harness fails ("mutation
smoke"). E4's implementer did this for golden drift; nobody did it for the
comparator until the reviewer.

## 10. Agents throw away their best evidence; make keeping it the default

Seen: E3 built a 110k-comparison R fuzz and gitignored it; E0 ran ~280 R
probes and kept only `# R:` comments. The verification happened; its
record didn't survive.
Try: a standing rule in every brief: any script that runs the reference
gets committed under `evidence/` (with versions), even if throwaway.
Cheap to keep, and it turns ad hoc probing into a re-runnable regression
check after an upstream change.

## 11. Planner stages should say who does them

Seen: E3's planner wrote "S1–S4 write the spec" as stages, so the
implementer (Sonnet) wrote the spec, not the planner (Opus).
Try: plan briefs say explicitly which stages the planner owns. If a
design artifact matters (spec, case schema), the planner writes it.

## 12. Keep a standing "port defect checklist" for reviewers

Seen: aliasing/mutability (E0b, E2), missing-value subscripts (E3),
length/recycling validation (E0), metadata dtype truncation on concat (E2).
None is an R-fidelity issue; all are found by probing Python contracts.
Try: a short checklist the reviewer always runs: mutate every accessor
result; feed NA/masked values to every index-taking method; wrong-length
vector args; mixed-dtype mcols through c(); empty inputs everywhere. Turn
it into a reusable probe script per class after one more round.

## 13. Fuzz is the step every approach converges on — schedule it

Seen: E0 (ad hoc probes), E2, E3 (fuzz) all fuzzed against R without
being told to; E2's and E3's fuzz found the same hotspots (seqlength
clipping, p* strand rules, nearest ties). Hand-authored cases (900 in E2)
missed them.
Try: make "random differential run against the reference" an explicit,
committed stage in every pipeline that has the reference available,
after hand-authored cases and before review. Commit seeds and failures as
cases.

## 14. A different-model adversary after "green" is the best-value stage

Seen: E1's red team (Sol, 15 min) found 52 divergences after the
implementer's 3,601 differential cases passed. Same-author cases share
the author's blind spots, even when generated at scale.
Try: in any pipeline with a reference, add a time-boxed adversary from a
different model family *after* the suite is green, with authority over
cases only (never src). Feed its failures to revise as fix-or-ADR.

## 15. "R is normative" needs an escape hatch for accidents

Seen: E1 reverse-engineered and cloned R's NA-weight coverage behavior,
which depends on which of two internal engines R picks. Every agent
treats "R is normative" literally.
Try: brief agents that when matching R requires reproducing a
path-dependent or engine-dependent quirk, they stop and write a
"contract or accident?" ADR stub instead of cloning. Batch the stubs for
a human (this is exactly where Hervé's judgment is needed).

## 16. Blind second adjudicator with frozen labels

Seen: Sol recorded its own class for each sampled cluster, hashed the
file, and only then read Opus's rationale. 92% agreement, and the 8%
disagreements were the most valuable output (the BiocPy over-excusal).
Try: any LLM-judged classification at scale gets a blind second judge
from a different model family on a stratified sample, labels frozen
before comparison. Report agreement per class; investigate classes that
shrink a denominator first.

## 17. Never fix a measurement by excluding cases; fix the instrument

Seen: "adapter_fault" was excluded from adjusted fidelity, which hid 148
real port failures behind a projection bug.
Try: when the instrument is wrong, repair it and re-measure. Exclusion
classes are for genuinely uncomparable cases only, and every exclusion
class gets the blind check (idea 16).

## 18. Invest in detection, not propagation

Seen: given a release note, all six approaches propagated an upstream
change correctly in 1–4 minutes. Only the approaches with re-runnable R
goldens detected it without being told.
Try: for any port of an evolving reference, the one non-negotiable artifact
is a cheap drift check run on every upstream upgrade (CI on Bioc devel).
Key it on *behavior* (re-run R), not on version strings (E1's cache missed
a same-version change).
