// Maintenance replay, propagation phase: an upstream behavior change lands; each approach
// propagates it with its own documented process. Same model (Sonnet) for all, in parallel.
// args: { wtRoot, rBin, news, ports: [{ exp, rAccess }] }

const IMPLEMENTER = "claude-bridge/claude-sonnet-5-5";

const brief = function (p) {
  const wt = args.wtRoot + "/replay-" + p.exp;
  const dir = "experiments/" + p.exp;
  const rLine = p.rAccess
    ? ["R/Bioconductor on this machine has been UPGRADED. The upgraded R is", args.rBin + "/Rscript.",
       "Start EVERY shell command with: export PATH=" + args.rBin + ":$PATH",
       "(so `Rscript`, `just` recipes and your harness all use the upgraded R). Never use",
       "/usr/local/bin/Rscript directly: that is the old installation."]
    : ["This project has no access to R (do not run R or Rscript); work from the release note."];
  return [
    "Maintenance task for the Python port in " + wt + "/" + dir + " (branch replay/" + p.exp + ").",
    "",
    "Upstream released new versions of IRanges/GenomicRanges. Release note:",
    args.news,
    "",
  ].concat(rLine, [
    "",
    "Bring the port in line with the new upstream behavior, following THIS PROJECT'S OWN",
    "documented maintenance process (read its README 'maintaining' section first, then NOTES,",
    "and any spec/ADRs/cases/goldens/harness it has). Update every artifact the process says to",
    "update (code, tests, specs, ADRs, cases, goldens, caches). Keep `uv run pytest` green.",
    "Write only inside " + dir + "/. Do not read other worktrees, branches or background/.",
    "Commit in stages (never squash/amend/rebase), messages naming the upstream change; end with",
    "Co-Authored-By: Claude <noreply@anthropic.com>. Foreground only. Autonomous. Terse.",
    "Reply with: what the process told you to do, what you changed (files), how you verified it",
    "against R (or why you could not), test counts, commits.",
  ]).join("\n");
};

const t0 = Date.now();
const results = await runs.all(args.ports.map(function (p) {
  return {
    key: "replay-" + p.exp,
    label: "Propagate upstream change " + p.exp,
    agent: "worker",
    model: IMPLEMENTER,
    cwd: args.wtRoot + "/replay-" + p.exp,
    timeoutMs: 90 * 60 * 1000,
    task: brief(p),
  };
}));

return {
  minutes_total: Math.round((Date.now() - t0) / 60000),
  results: args.ports.map(function (p, i) {
    return { exp: p.exp, runId: results[i].runId, ok: results[i].ok, output: results[i].output };
  }),
};
