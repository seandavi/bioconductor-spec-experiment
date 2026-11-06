// One experiment: plan (Opus) -> implement (Sonnet) -> [redteam (Sol)] -> review (Sol) -> revise (Sonnet).
// Same pipeline and models for every experiment; only the approach brief and flags differ.
//
// args: {
//   exp: "e1-differential",
//   worktree: "/abs/path/to/worktree",
//   approach: "<approach brief>",
//   rAccess: true,    // false: agents may not run R or read R sources/docs
//   redteam: false,   // true: adversarial differential-case stage before review
// }
//
// History: E0 ran v1 (commit d72e8bd): reviewer was the read-only `reviewer`
// agent (no shell), no terseness rule. v2: reviewer is a no-edit `worker`
// with shell, terse-writing rule, rAccess and redteam flags.

const PLANNER = "claude-bridge/claude-opus-5-5";
const IMPLEMENTER = "claude-bridge/claude-sonnet-5-5";
const REVIEWER = "openai-codex/gpt-6.1-sol";

const dir = "experiments/" + args.exp;
const pkg = "ranges_" + args.exp.split("-")[0];
const rAccess = args.rAccess !== false;

const rRule = rAccess
  ? [
      "- R 4.6 with Bioc 3.24 (IRanges, GenomicRanges, testthat) is installed; `Rscript`",
      "  is on PATH. R is the normative reference. How you may use R is set by your",
      "  approach below.",
    ]
  : [
      "- R is the normative reference, but in this experiment you have NO access to it:",
      "  do not run `R` or `Rscript`, and do not read R package sources, help files or",
      "  vignettes on this machine. Work from your own knowledge of IRanges/GenomicRanges.",
    ];

const common = [
  "You are part of a controlled experiment comparing processes for porting",
  "Bioconductor's IRanges/GRanges to Python. Your experiment is " + args.exp + ".",
  "",
  "Repository: git worktree at " + args.worktree + " on branch exp/" + args.exp + ".",
  "Read TARGET.md (scope, done criteria) and PROTOCOL.md (rules) at the worktree root first.",
  "",
  "Rules that apply to you:",
  "- Write only inside " + dir + "/ in this worktree. Python package: uv project,",
  "  src layout, import name `" + pkg + "`. Type hints throughout.",
  "- Do NOT read other branches, other worktrees (anything under",
  "  bioconductor-spec-experiment.worktrees/ except your own), the main checkout,",
  "  or any background/ directory. Do not run `git log --all` or inspect other refs.",
  "- Do NOT install, import, vendor or read BiocPy (`iranges`, `genomicranges`,",
  "  `biocframe`, `biocutils`, `dolomite-*` on PyPI) or its source. numpy/polars etc. are fine.",
].concat(rRule, [
  "- Commit at every meaningful stage with `git commit` (never squash, never amend,",
  "  no rebases). Messages say what changed and why, especially redesigns.",
  "  End each message with: Co-Authored-By: Claude <noreply@anthropic.com>",
  "- Work in the foreground: run tests and commands directly and wait for them.",
  "  Never start background jobs and wait for notifications.",
  "- Fully autonomous: do not ask the human anything. Record judgment calls in",
  "  " + dir + "/NOTES.md (append-only, timestamped).",
  "- Write tersely everywhere: docs, specs, ADRs, NOTES, commit messages, replies.",
  "  Short sentences, no filler, no restating the obvious; prefer bullets and tables.",
  "  Anyone should be able to skim any single document in two minutes.",
  "",
  "Approach for this experiment:",
  args.approach,
]).join("\n");

const t0 = Date.now();
const times = {};
let tStage = Date.now();
const lap = function (name) { times[name] = Date.now() - tStage; tStage = Date.now(); };

const plan = await runs.run("plan", {
  label: "Plan " + args.exp,
  agent: "worker",
  model: PLANNER,
  cwd: args.worktree,
  timeoutMs: 45 * 60 * 1000,
  task: common + "\n\nYOUR STAGE: PLAN. You are the planner; a different agent implements.\n" +
    "Write " + dir + "/PLAN.md: architecture (data model, core algorithms, module layout),\n" +
    "the artifacts your approach calls for and how they are maintained, a staged work\n" +
    "breakdown that an implementer can follow commit by commit, and the main risks.\n" +
    "Do not implement the package. You may scaffold the uv project (pyproject, empty src\n" +
    "package, tests dir) so the implementer starts clean. Commit PLAN.md and the scaffold.\n" +
    "Reply with a short summary of the plan and the commit hashes.",
});
lap("plan");

let impl = await runs.run("implement", {
  label: "Implement " + args.exp,
  agent: "worker",
  model: IMPLEMENTER,
  cwd: args.worktree,
  timeoutMs: 3 * 60 * 60 * 1000,
  task: common + "\n\nYOUR STAGE: IMPLEMENT. Follow " + dir + "/PLAN.md stage by stage,\n" +
    "committing after each stage. You may deviate from the plan when you find it wrong;\n" +
    "record why in NOTES.md and the commit message. Stop when the TARGET.md done criteria\n" +
    "are met (all required operations public, `uv run pytest` green in " + dir + ",\n" +
    "README.md written, everything committed). Reply with: what is done, what is not,\n" +
    "final test counts, and the list of commits.",
});
lap("implement");

let redteam = null;
if (args.redteam) {
  redteam = await runs.run("redteam", {
    label: "Red-team " + args.exp,
    agent: "worker",
    model: REVIEWER,
    cwd: args.worktree,
    timeoutMs: 90 * 60 * 1000,
    task: common + "\n\nYOUR STAGE: RED TEAM. You did not write this port; your job is to\n" +
      "break it. Using the experiment's differential harness (extend it if needed), derive\n" +
      "test cases most likely to expose behavior where the Python port differs from R:\n" +
      "edge cases, argument combinations, strand/seqinfo interactions, empty and zero-width\n" +
      "inputs, ordering, names and mcols propagation. Run each in R and Python.\n" +
      "Commit the cases (all of them, passing and failing) in stages. Do NOT modify\n" +
      "src/ or the existing tests' expectations; mark failing cases as expected failures\n" +
      "with a reason so the suite stays green. Reply with a table: case id, what R does,\n" +
      "what Python does, and your view of whether it is a bug, an acceptable difference,\n" +
      "or a judgment call needing an ADR.",
  });
  lap("redteam");
}

const review = await runs.run("review", {
  label: "Review " + args.exp,
  agent: "worker",
  model: REVIEWER,
  cwd: args.worktree,
  timeoutMs: 60 * 60 * 1000,
  output: "reviews/" + args.exp + "-review.md",
  task: common + "\n\nYOUR STAGE: REVIEW. You may run commands (tests, the experiment's\n" +
    "harnesses, small probe scripts under /tmp) but must NOT edit, create or commit\n" +
    "anything in the worktree. Review " + dir + "/ at HEAD of exp/" + args.exp + ".\n" +
    "Scope your review the same way for every experiment:\n" +
    "1. Completeness against TARGET.md required scope (list missing or partial operations).\n" +
    "2. Correctness against the experiment's OWN stated behavior: its tests, docs, and any\n" +
    "   specs/ADRs/cases it wrote. Run `uv run pytest` and any harness it ships.\n" +
    "3. Whether the experiment actually followed its stated approach.\n" +
    "4. Code quality, design and maintainability.\n" +
    "Do NOT do your own independent R-vs-Python fidelity hunt (that is evaluated\n" +
    "separately, after all experiments). Only raise R behavior where the experiment's own\n" +
    "artifacts make a claim about it that is internally inconsistent.\n" +
    "Return numbered findings, each with severity (blocker/major/minor), file:line, and\n" +
    "a concrete fix.",
});
lap("review");

let reviseInput = "REVIEW:\n" + review.output;
if (redteam) {
  reviseInput = "RED-TEAM REPORT (differential cases now committed; failing ones marked " +
    "expected-fail):\n" + redteam.output + "\n\n" + reviseInput;
}

impl = await runs.run("revise", {
  label: "Revise " + args.exp,
  resume: impl.runId,
  timeoutMs: 2 * 60 * 60 * 1000,
  task: "Independent agents examined your work (reports below). Address every blocker\n" +
    "and major finding, and minor ones where cheap." +
    (redteam ? " For each failing red-team case: fix the port to match R and drop the\n" +
      "expected-fail mark, or, if matching R is wrong or not worth it, write a short ADR\n" +
      "and keep the case documenting the divergence." : "") +
    "\nFor any finding you reject, say why in NOTES.md. Commit in stages (one commit per\n" +
    "finding or tight group, messages referencing finding/case ids). Keep `uv run pytest`\n" +
    "green. Reply with what you fixed, what you rejected (and why), and the commits.\n\n" +
    reviseInput,
});
lap("revise");
times.total = Date.now() - t0;

return {
  exp: args.exp,
  times_ms: times,
  plan: { runId: plan.runId, ok: plan.ok, output: plan.output },
  redteam: redteam ? { runId: redteam.runId, ok: redteam.ok, output: redteam.output } : null,
  review: { runId: review.runId, ok: review.ok, outputReference: review.outputReference },
  final: { runId: impl.runId, ok: impl.ok, output: impl.output },
};
