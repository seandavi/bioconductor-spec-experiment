// Eval phase 2: one adapter per port, each in its own worktree.
// write (Sonnet) x7 in parallel -> thinness review (Sol) x7 -> revise (Sonnet, resumed) x7.
// args: { wtRoot: "/abs/worktrees", ports: [{ port, dir, api }] }
//   dir: the port's uv project dir (absolute); api: one line on where its public API is documented.

const WRITER = "claude-bridge/claude-sonnet-5-5";
const REVIEWER = "openai-codex/gpt-6.1-sol";

const brief = function (p) {
  const wt = args.wtRoot + "/eval-adapter-" + p.port;
  return [
    "Write the eval adapter for port `" + p.port + "`.",
    "",
    "Worktree: " + wt + " (branch eval-adapter/" + p.port + "). Work from " + wt + "/eval.",
    "Read eval/ADAPTERS.md first and follow it exactly; also spec/PROJECTION.md, spec/ops.json,",
    "spec/CASES.md. You write ONE file: eval/adapters/" + p.port + ".py (plus a short",
    "eval/adapters/" + p.port + ".md noting unsupported ops and convention translations).",
    "",
    "The port: " + p.dir + " (frozen, READ-ONLY). " + p.api,
    "Read its README, docstrings and public names in src/. Do NOT read its tests, cases,",
    "goldens or ADRs; do NOT read other adapters, other experiments or background/.",
    "",
    "THIN means: translate R calling conventions to the port's API and project results to",
    "PROJECTION.md (index bases, NA representation, argument names, object construction).",
    "Never add semantic logic that changes what the port computes: no fixing, sorting,",
    "deduplicating, clipping or re-deriving results. If the port computes something wrong,",
    "the score must show it. If the port lacks an op or argument, declare it unsupported.",
    "",
    "Dev loop: `just roundtrip " + p.port + "` and `just adapter-check " + p.port + "` only.",
    "Do NOT run `just run`, `just report`, or read build/expected.jsonl or any results:",
    "you must not see R's expected values or scores.",
    "Done when adapter-check shows no adapter_error you can remove without semantic logic",
    "and roundtrip passes for every input kind the port supports.",
    "",
    "Commit in stages on eval-adapter/" + p.port + " (never squash/amend/rebase); end messages",
    "with Co-Authored-By: Claude <noreply@anthropic.com>. Foreground only. Autonomous.",
    "Terse prose. Reply with: ops supported/unsupported counts, remaining adapter_errors and",
    "why, conventions translated, commits.",
  ].join("\n");
};

const writes = await runs.all(args.ports.map(function (p) {
  return {
    key: "write-" + p.port,
    label: "Write adapter " + p.port,
    agent: "worker",
    model: WRITER,
    cwd: args.wtRoot + "/eval-adapter-" + p.port,
    timeoutMs: 2 * 60 * 60 * 1000,
    task: brief(p),
  };
}));

const reviews = await runs.all(args.ports.map(function (p) {
  return {
    key: "review-" + p.port,
    label: "Review adapter " + p.port,
    agent: "worker",
    model: REVIEWER,
    cwd: args.wtRoot + "/eval-adapter-" + p.port,
    timeoutMs: 45 * 60 * 1000,
    output: "reviews/eval-adapter-" + p.port + "-review.md",
    task: "Review the eval adapter eval/adapters/" + p.port + ".py in worktree " +
      args.wtRoot + "/eval-adapter-" + p.port + " against eval/ADAPTERS.md and the port's public\n" +
      "API (" + p.dir + ", read-only). Do NOT edit or commit anything in the worktree; /tmp probes\n" +
      "are fine. Do NOT run `just run`/`just report` or read build/expected.jsonl.\n" +
      "Check, in order:\n" +
      "1. THINNESS (most important): any logic that changes what the port computes (fixing,\n" +
      "   sorting, deduplicating, clipping, re-deriving, catching port errors and substituting\n" +
      "   values) would inflate this port's score. List every instance.\n" +
      "2. Translation correctness: index bases, NA, argument names/defaults, strand/seqinfo\n" +
      "   construction, projection shape. A translation bug would deflate the score unfairly.\n" +
      "3. Unsupported declarations: is anything declared unsupported that the port's public API\n" +
      "   actually provides (or vice versa)?\n" +
      "Run `just roundtrip " + p.port + "` and `just adapter-check " + p.port + "`.\n" +
      "Return numbered findings with severity, file:line and a concrete fix. Terse.",
  };
}));

const revises = await runs.all(args.ports.map(function (p, i) {
  return {
    key: "revise-" + p.port,
    label: "Revise adapter " + p.port,
    resume: writes[i].runId,
    timeoutMs: 60 * 60 * 1000,
    task: "An independent reviewer examined your adapter. Fix every thinness finding and every\n" +
      "translation bug; reject only with evidence (note it in eval/adapters/" + p.port + ".md).\n" +
      "Same rules as before (no `just run`, no expected values). Commit per finding. Reply with\n" +
      "fixes, rejections, final adapter-check summary, commits.\n\nREVIEW:\n" + reviews[i].output,
  };
}));

return args.ports.map(function (p, i) {
  return {
    port: p.port,
    write: { runId: writes[i].runId, ok: writes[i].ok },
    review: { runId: reviews[i].runId, ok: reviews[i].ok, outputReference: reviews[i].outputReference },
    revise: { runId: revises[i].runId, ok: revises[i].ok, output: revises[i].output },
  };
});
