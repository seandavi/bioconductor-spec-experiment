// Eval phase 1: build the cross-port fidelity harness.
// plan (Opus) -> build (Sonnet) -> harness review (Sol) -> revise (Sonnet).
// args: { worktree: "/abs/eval/worktree", wtRoot: "/abs/worktrees dir", ports: "<port table text>" }

const PLANNER = "claude-bridge/claude-opus-5-5";
const BUILDER = "claude-bridge/claude-sonnet-5-5";
const REVIEWER = "openai-codex/gpt-6.1-sol";

const common = [
  "You are building the EVALUATION harness for a controlled experiment: six independent",
  "Python ports of Bioconductor IRanges/GRanges (plus BiocPy as a bonus) will be scored on",
  "fidelity to R. R is normative.",
  "",
  "Repository: git worktree at " + args.worktree + " on branch `eval`. Read TARGET.md and",
  "PROTOCOL.md at its root. Write only inside eval/ in this worktree.",
  "",
  "Ports (frozen; READ-ONLY: never edit, commit, or run uv lock/sync that changes them):",
  args.ports,
  "",
  "Hard requirements:",
  "- Expected values come from running R LIVE through the eval's own R runner. Never use any",
  "  experiment's goldens, caches or expected values. Caching the eval runner's own R output",
  "  (keyed by R + package versions and case hash) is fine.",
  "- Case pool, every case tagged with `source`:",
  "  (a) converted from the E1, E2 and E4 suites: inputs and calls only, expected values",
  "      dropped. E1: " + args.wtRoot + "/e1-differential/experiments/e1-differential/tests/differential/",
  "      (YAML in R vocabulary; red-team cases tagged `e1-redteam` where identifiable).",
  "      E2: " + args.wtRoot + "/e2-conformance/experiments/e2-conformance/conformance/cases/.",
  "      E4: " + args.wtRoot + "/e4-primitive-cases/experiments/e4-primitive-cases/cases/.",
  "      Cases that cannot be expressed in the eval schema are counted and listed, not hidden.",
  "  (b) `heldout`: random cases from a seeded generator written from TARGET.md scope only,",
  "      covering every required op family (strands, seqinfo, zero-width, empty inputs,",
  "      names, mcols), with seeds no experiment has used. Write it without looking at any",
  "      port's tests or cases.",
  "- One normalized result projection in R conventions (1-based indices, NA as null, typed",
  "  vectors) per result kind (IRanges, GRanges, Hits, integer/logical vectors, Rle/coverage,",
  "  Seqinfo, error).",
  "- Comparator: exact by default. Allowed relaxations, each documented in one place:",
  "  Hits compared as a set only where R's order is not a contract (justify);",
  "  select=\"arbitrary\" passes if the port's hit is one of R's select=\"all\" hits;",
  "  both-error = pass (messages not compared). R error vs port success, or the reverse,",
  "  is a mismatch.",
  "- Outcome categories per (case, port): pass, mismatch, adapter_error (the adapter itself",
  "  failed: TypeError in translation, unknown API), unsupported (adapter declares the op",
  "  absent from the port), port_error_unexpected (port raised where R did not).",
  "- Port runner contract: each port runs in ITS OWN venv via",
  "  `uv run --project <port dir> python <eval runner> --port <name> ...`, so the runner and",
  "  adapter base code may import ONLY the Python stdlib (+ numpy, present in every port) and",
  "  the port. JSON/JSONL on disk between stages. BiocPy gets its own small uv project under",
  "  eval/ports/biocpy pinned to iranges==0.7.3, genomicranges==0.9.0.",
  "- Adapters are NOT written in this phase except two self-tests: `r_echo` (returns R's own",
  "  projection, must score 100%) and `mutant` (r_echo with a deterministic perturbation, must",
  "  fail the perturbed cases). Define the adapter interface clearly in eval/ADAPTERS.md so",
  "  seven fresh agents can each write one adapter from it plus one port's public API.",
  "- Report: matrix of port x op family x source with counts per outcome; per-case JSONL;",
  "  a `justfile` with `just oracle`, `just run PORT`, `just report`.",
  "- Commit in stages; never squash/amend/rebase. End messages with",
  "  Co-Authored-By: Claude <noreply@anthropic.com>. Work in the foreground. Fully",
  "  autonomous. Write tersely (bullets/tables; any doc skimmable in two minutes).",
].join("\n");

const times = {};
let t = Date.now();
const lap = function (k) { times[k] = Date.now() - t; t = Date.now(); };

const plan = await runs.run("plan", {
  label: "Plan eval harness",
  agent: "worker",
  model: PLANNER,
  cwd: args.worktree,
  timeoutMs: 45 * 60 * 1000,
  task: common + "\n\nYOUR STAGE: PLAN + DESIGN ARTIFACTS. You own the design. Write\n" +
    "eval/EVAL-PLAN.md (architecture, staged build breakdown, risks) AND write yourself the\n" +
    "artifacts that define the contract: the case schema, the projection spec per result kind,\n" +
    "the comparator rules, and eval/ADAPTERS.md. Inspect the three source suites' formats to\n" +
    "design the converters. Do not build the runners. Commit. Reply with a short summary.",
});
lap("plan");

let build = await runs.run("build", {
  label: "Build eval harness",
  agent: "worker",
  model: BUILDER,
  cwd: args.worktree,
  timeoutMs: 3 * 60 * 60 * 1000,
  task: common + "\n\nYOUR STAGE: BUILD. Implement eval/EVAL-PLAN.md stage by stage against\n" +
    "the planner's schema/projection/comparator/ADAPTERS.md (change them only with a NOTES\n" +
    "entry saying why). Done when: converters run over all three suites (with counts of\n" +
    "converted/unconvertible per suite), the held-out generator emits cases for every op\n" +
    "family, `just oracle` fills R outputs for the whole pool, `r_echo` scores 100% and\n" +
    "`mutant` fails exactly its perturbed cases, `just report` renders, and harness unit tests\n" +
    "pass. Reply with counts per source, oracle runtime, and commits.",
});
lap("build");

const review = await runs.run("review", {
  label: "Review eval harness",
  agent: "worker",
  model: REVIEWER,
  cwd: args.worktree,
  timeoutMs: 60 * 60 * 1000,
  output: "reviews/eval-harness-review.md",
  task: common + "\n\nYOUR STAGE: REVIEW THE ORACLE (no edits or commits in the worktree;\n" +
    "probes under /tmp are fine). This harness will referee six ports, so its defects become\n" +
    "every port's score. Check: converter fidelity (sample cases from each suite; are inputs\n" +
    "and calls faithful?), projection completeness per result kind, comparator strictness\n" +
    "(perturb values, names, order, NA, types and confirm failure; confirm the relaxations\n" +
    "are justified and narrow), held-out coverage of every TARGET op family, the outcome\n" +
    "categorisation (can an adapter bug masquerade as a port mismatch?), and whether\n" +
    "ADAPTERS.md is sufficient for a fresh agent. Return numbered findings with severity,\n" +
    "file:line and a concrete fix.",
});
lap("review");

build = await runs.run("revise", {
  label: "Revise eval harness",
  resume: build.runId,
  timeoutMs: 90 * 60 * 1000,
  task: "An independent reviewer examined the harness. Fix every blocker and major finding,\n" +
    "minor ones where cheap; reject only with evidence, recorded in eval/NOTES.md. Commit per\n" +
    "finding. Re-run r_echo (100%) and mutant (fails exactly the perturbed cases) at the end.\n" +
    "Reply with fixes, rejections, and commits.\n\nREVIEW:\n" + review.output,
});
lap("revise");

return {
  times_ms: times,
  plan: { runId: plan.runId, output: plan.output },
  review: { runId: review.runId, outputReference: review.outputReference },
  final: { runId: build.runId, ok: build.ok, output: build.output },
};
