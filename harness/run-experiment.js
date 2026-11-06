// One experiment: plan (Opus) -> implement (Sonnet) -> review (Sol) -> revise (Sonnet).
// Same pipeline and models for every experiment; only the approach brief differs.
//
// args: { exp: "e0-control", worktree: "/abs/path/to/worktree", approach: "<approach brief>" }

const PLANNER = "claude-bridge/claude-opus-5-5";
const IMPLEMENTER = "claude-bridge/claude-sonnet-5-5";
const REVIEWER = "openai-codex/gpt-6.1-sol";

const dir = "experiments/" + args.exp;
const pkg = "ranges_" + args.exp.split("-")[0];

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
  "  `biocframe`, `biocutils` on PyPI) or its source. numpy/polars etc. are fine.",
  "- R 4.6 with Bioc 3.24 (IRanges, GenomicRanges, testthat) is installed; `Rscript`",
  "  is on PATH. R is the normative reference. How you may use R is set by your",
  "  approach below.",
  "- Commit at every meaningful stage with `git commit` (never squash, never amend",
  "  pushed work, no rebases). Messages say what changed and why, especially",
  "  redesigns. End each message with: Co-Authored-By: Claude <noreply@anthropic.com>",
  "- Work in the foreground: run tests and commands directly and wait for them.",
  "  Never start background jobs and wait for notifications.",
  "- Fully autonomous: do not ask the human anything. Record judgment calls in",
  "  " + dir + "/NOTES.md (append-only, timestamped).",
  "",
  "Approach for this experiment:",
  args.approach,
].join("\n");

const t0 = Date.now();
const times = {};

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
times.plan = Date.now() - t0;

const t1 = Date.now();
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
times.implement = Date.now() - t1;

const t2 = Date.now();
const review = await runs.run("review", {
  label: "Review " + args.exp,
  agent: "reviewer",
  model: REVIEWER,
  cwd: args.worktree,
  timeoutMs: 60 * 60 * 1000,
  output: "reviews/" + args.exp + "-review.md",
  task: common + "\n\nYOUR STAGE: REVIEW (read-only; do not edit or commit).\n" +
    "Review " + dir + "/ at HEAD of exp/" + args.exp + ". Scope your review the same way\n" +
    "for every experiment:\n" +
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
times.review = Date.now() - t2;

const t3 = Date.now();
impl = await runs.run("revise", {
  label: "Revise " + args.exp,
  resume: impl.runId,
  timeoutMs: 90 * 60 * 1000,
  task: "An independent reviewer reviewed your work. Address every blocker and major\n" +
    "finding, and minor ones where cheap. For any finding you reject, say why in\n" +
    "NOTES.md. Commit in stages as before (one commit per finding or tight group,\n" +
    "messages referencing the finding numbers). Keep `uv run pytest` green.\n" +
    "Reply with which findings you fixed, rejected (and why), and the commits.\n\n" +
    "REVIEW:\n" + review.output,
});
times.revise = Date.now() - t3;
times.total = Date.now() - t0;

return {
  exp: args.exp,
  times_ms: times,
  plan: { runId: plan.runId, ok: plan.ok, output: plan.output },
  review: { runId: review.runId, ok: review.ok, outputReference: review.outputReference },
  final: { runId: impl.runId, ok: impl.ok, output: impl.output },
};
