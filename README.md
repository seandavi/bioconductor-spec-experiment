# bioconductor-spec-experiment

How much process does a faithful Python port of a Bioconductor class
need? Five independent agentic ports of IRanges/GRanges, each using a
different approach (plain port, differential testing, formal
conformance suite, spec + ADRs, primitive shared cases), compared on
design, velocity, fidelity to R and maintenance cost.

- [TARGET.md](TARGET.md) — what gets built, and when it counts as done
- [PROTOCOL.md](PROTOCOL.md) — the experiments, isolation rules, metrics
- `background/` — the Bioconductor Zulip discussion that prompted this
