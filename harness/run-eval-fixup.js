// Eval phase 4: fix adapter/harness faults found in adjudication, re-score, re-adjudicate, re-check.
// fix (Sonnet) -> re-adjudicate (Opus, resumed) -> verify (Sol).
// args: { worktree, adjudicatorRunId, spotcheckPath }

const FIXER = "claude-bridge/claude-sonnet-5-5";
const CHECKER = "openai-codex/gpt-6.1-sol";

const fix = await runs.run("fix", {
  label: "Fix adapter and harness faults",
  agent: "worker",
  model: FIXER,
  cwd: args.worktree,
  timeoutMs: 2 * 60 * 60 * 1000,
  task: "Eval worktree " + args.worktree + " (branch eval). Read eval/adjudication/ADJUDICATION.md\n" +
    "(sections on adapter_fault and harness_fault) and the independent spot-check at\n" +
    args.spotcheckPath + " (read it fully).\n" +
    "Fix, one commit each, in eval/adapters/ and eval/runner/ only:\n" +
    "1. Every adapter_fault and the harness_fault listed in ADJUDICATION.md, with the stated fix.\n" +
    "2. Spot-check items: S031 (biocpy select= projection: representation fix only; never\n" +
    "   recompute hits); S050 (translate R incomparables=FALSE to the native generic's\n" +
    "   equivalent); S051/S056/S093/S104 (vector forms the native public API does not take\n" +
    "   -> Unsupported, not forwarded); S014/S029 (zero-column mcols: project with the owner's\n" +
    "   length/names, a representation fix).\n" +
    "Adapters stay THIN: never change what a port computes. Then run, from eval/:\n" +
    "`just selftest` (r_echo 100%, mutant 0 disagreements), `just run <port>` for e0 e0b e1 e2 e3\n" +
    "e4 biocpy, `just report`. Commit the new results. Ports are frozen and read-only.\n" +
    "Commits: never squash/amend; end with Co-Authored-By: Claude <noreply@anthropic.com>.\n" +
    "Foreground only. Autonomous. Terse. Reply with fixes+commits and the new headline table.",
});

const readj = await runs.run("readjudicate", {
  label: "Re-adjudicate after fixes",
  resume: args.adjudicatorRunId,
  timeoutMs: 2 * 60 * 60 * 1000,
  task: "The adapter/harness faults you listed have been fixed and every port re-scored (see\n" +
    "eval/results/ and recent commits on `eval`). An independent spot-check of your work is at\n" +
    args.spotcheckPath + "; read it fully. Now:\n" +
    "1. Re-run cluster.py --write and build.py --check on the new results; add or amend rules in\n" +
    "   verdicts.py for every unmatched or changed signature.\n" +
    "2. Apply the spot-check's split requests: S031 residual biocpy select= failures (port_bug,\n" +
    "   not adapter_fault), S041/S051/S056/S050 mixed clusters split by cause, S073 rationale\n" +
    "   rewritten, S116 port_bug. Where a port both hits an R accident AND has its own defect on\n" +
    "   the same case, classify port_bug (do not let r_accident excuse a separate port defect).\n" +
    "3. adapter_fault should now be ~0; any that remain, list with cause.\n" +
    "4. Rewrite ADJUDICATION.md with the new counts, adjusted held-out fidelity, and a short\n" +
    "   'changes since first adjudication' section. Commit in stages. Reply with the new per-port\n" +
    "   class table, adjusted fidelity, and commits.",
});

const verify = await runs.run("verify", {
  label: "Verify re-adjudication",
  agent: "worker",
  model: CHECKER,
  cwd: args.worktree,
  timeoutMs: 60 * 60 * 1000,
  output: "reviews/eval-readjudication-verify.md",
  task: "Eval worktree " + args.worktree + ". You did the earlier spot-check (" + args.spotcheckPath +
    "). The adapter/harness faults have been fixed, ports re-scored, and the adjudication\n" +
    "revised (eval/adjudication/, recent commits). Do not edit or commit; /tmp probes fine.\n" +
    "Check: (1) every disagreement and split request from your spot-check is resolved or\n" +
    "explicitly answered; (2) the adapter fixes are thin (no semantic logic) and raised no\n" +
    "port's score by changing what it computes; (3) r_accident is not excusing separate port\n" +
    "defects; (4) re-run 2 cases from each of the 5 largest changed clusters. Report a verdict\n" +
    "(accept / accept with notes / reject) and any remaining issues. Terse.",
});

return {
  fix: { runId: fix.runId, ok: fix.ok, output: fix.output },
  readjudicate: { runId: readj.runId, ok: readj.ok, output: readj.output },
  verify: { runId: verify.runId, ok: verify.ok, outputReference: verify.outputReference },
};
