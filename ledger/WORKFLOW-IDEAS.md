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
