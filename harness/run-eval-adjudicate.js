// Eval phase 3: adjudicate mismatch clusters. Opus classifies; Sol independently re-classifies a sample.
// args: { worktree: "/abs/eval worktree", wtRoot: "/abs/worktrees" }

const ADJUDICATOR = "claude-bridge/claude-opus-5-5";
const CHECKER = "openai-codex/gpt-6.1-sol";

const ctx = [
  "Context: six independent Python ports of Bioconductor IRanges/GRanges (e0, e0b, e1, e2, e3, e4)",
  "plus BiocPy (iranges 0.7.3, genomicranges 0.9.0) were scored against live R on 9,650 cases.",
  "Eval worktree: " + args.worktree + " (branch eval). Read eval/README.md, eval/spec/COMPARATOR.md,",
  "eval/spec/PROJECTION.md, eval/results/report.md. Scored per-case records: eval/results/scored/.",
  "Cases: eval/build/resolved.jsonl; R expected: eval/build/expected.jsonl; adapters: eval/adapters/.",
  "Ports (READ-ONLY, frozen; you may read anything in them now, including ADRs, NOTES, README):",
  "  e0 " + args.wtRoot + "/e0-control/experiments/e0-control",
  "  e0b " + args.wtRoot + "/e0b-no-r/experiments/e0b-no-r",
  "  e1 " + args.wtRoot + "/e1-differential/experiments/e1-differential",
  "  e2 " + args.wtRoot + "/e2-conformance/experiments/e2-conformance",
  "  e3 " + args.wtRoot + "/e3-spec-adr/experiments/e3-spec-adr",
  "  e4 " + args.wtRoot + "/e4-primitive-cases/experiments/e4-primitive-cases",
  "  biocpy: eval/ports/biocpy venv (site-packages iranges, genomicranges)",
  "Rscript (R 4.6, Bioc 3.24) is available; re-run any case in R or in a port (`uv run --project`).",
  "Never edit ports. Write tersely: tables and bullets.",
].join("\n");

const t0 = Date.now();

const adj = await runs.run("adjudicate", {
  label: "Adjudicate mismatch clusters",
  agent: "worker",
  model: ADJUDICATOR,
  cwd: args.worktree,
  timeoutMs: 3 * 60 * 60 * 1000,
  task: ctx + "\n\nYOUR TASK: adjudicate every non-pass, non-unsupported outcome (mismatch,\n" +
    "port_error_unexpected) for every port. Write only under eval/adjudication/.\n" +
    "1. Cluster by a reproducible signature (port, op, mismatch kind, decisive input feature);\n" +
    "   commit the clustering script. Aim for clusters a person can reason about (tens per port).\n" +
    "2. For each cluster, examine >=3 cases (all if <=3): R's result, the port's result, the\n" +
    "   adapter path, and the port's docs/ADRs. Classify as exactly one of:\n" +
    "   - port_bug: the port computes something R does not, with no declared reason.\n" +
    "   - declared_divergence: the port documents this difference (cite ADR/README/NOTES).\n" +
    "   - r_accident: R's behavior is path/engine-dependent, internally inconsistent, or\n" +
    "     defective, so matching it is arguably not a contract (give evidence, e.g. two\n" +
    "     equivalent R calls disagreeing). Be strict: 'surprising' is not an accident.\n" +
    "   - adapter_fault: the adapter mistranslated (cite adapter line); the port is not at fault.\n" +
    "   - harness_fault: converter/oracle/projection/comparator defect (cite file:line).\n" +
    "3. Write eval/adjudication/clusters.jsonl (signature, port, n, case ids, class, evidence,\n" +
    "   one-line rationale) and eval/adjudication/ADJUDICATION.md: per-port table of counts by\n" +
    "   class; adjusted held-out fidelity = pass / (supported - r_accident - declared -\n" +
    "   adapter_fault - harness_fault); the 10 most consequential port_bug clusters overall;\n" +
    "   every adapter_fault and harness_fault with the fix needed (do not apply fixes).\n" +
    "4. Specifically answer, with case evidence: how does each port (and BiocPy) order\n" +
    "   seqlevels when a GRanges is built without explicit seqinfo, vs R (order of first\n" +
    "   appearance)? Which ports implement circular-sequence semantics at all?\n" +
    "Commit in stages (never squash/amend); end messages with Co-Authored-By: Claude\n" +
    "<noreply@anthropic.com>. Foreground only. Autonomous. Reply with the per-port class table,\n" +
    "adjusted fidelity, adapter/harness faults found, and commits.",
});

const check = await runs.run("spotcheck", {
  label: "Spot-check adjudication",
  agent: "worker",
  model: CHECKER,
  cwd: args.worktree,
  timeoutMs: 90 * 60 * 1000,
  output: "reviews/eval-adjudication-spotcheck.md",
  task: ctx + "\n\nYOUR TASK: independently check the adjudication in eval/adjudication/ (do not\n" +
    "edit or commit anything; /tmp probes fine). Select: ALL biocpy clusters, plus for each of\n" +
    "e0..e4 the 3 largest clusters and 2 random clusters (seed 7; state which). For each,\n" +
    "re-examine >=2 cases yourself (run R and the port) BEFORE reading the adjudicator's\n" +
    "rationale, record your own class, then compare. Report: agreement rate overall and per\n" +
    "class, every disagreement with your evidence, and any cluster whose signature lumps\n" +
    "together cases that deserve different classes. Pay special attention to r_accident\n" +
    "(over-use would flatter ports) and adapter_fault (over-use would excuse ports).",
});

return {
  minutes: Math.round((Date.now() - t0) / 60000),
  adjudicate: { runId: adj.runId, ok: adj.ok, output: adj.output },
  spotcheck: { runId: check.runId, ok: check.ok, outputReference: check.outputReference },
};
