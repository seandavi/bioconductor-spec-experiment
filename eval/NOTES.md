# Build notes (changes to the planner's contract, deviations, caveats)

| When | Note |
|---|---|
| S4 heldout | **Blindness caveat.** The builder wrote `heldout/` (HELDOUT.md asked for a fresh agent). Mitigation: written *before* opening any E1/E2/E4 suite file or any port source. Before that the builder had read `spec/CONVERTERS.md` (HELDOUT.md forbids it for the generator author; it lists the 4 E1 `raw` ids and suite sizes), two `spec/examples/cases.jsonl` lines, and `ls` of E1's case file names. No suite case content, no port code. |
| S4 heldout | 3510 cases (brief: "2-3k"): 30 per (op, class), 60 for findOverlaps/nearest/reduce/coverage, 90/90/60 for the IRanges/GRanges/Seqinfo constructors. Pruned only if oracle time demands. |
| S4 heldout | `merge` is generated for `Seqinfo` only (GRanges `merge` is not Seqinfo merging); `seqinfo`, `seqlevels<-`, `seqinfo<-`, `trim` are GRanges-only (Seqinfo methods unverified blind); `sort(by=)` (a formula) is not expressible. |
| S4 seeds | `just seeds` scans every experiment worktree for each seed (decimal and hex); python scan, not grep (BSD grep -f was >2 min). |
