# Background

This experiment grew out of a discussion on the Bioconductor community Zulip
(November 2026, following a conference talk) about whether the core
Bioconductor classes (IRanges, GRanges, SummarizedExperiment and relatives)
should have formal, implementation-independent specifications, in the style
of W3C specs, so that R, Python and other implementations share one contract.

Positions raised in that thread, paraphrased:

- Write behavioral specs for the core classes, independent of any one
  implementation, so ports such as BiocPy have a reference to work from.
- Use AI to draft specs (JSON Schema, YAML or prose) from the existing R/S4
  code, with humans reviewing and maintaining them.
- Specs alone are not enough: conformance tests tied to decisions (ADRs) are
  what pin behavior down, and they can be derived from existing test suites.
- Beware overspecification: semantically equivalent results should not have
  to be identical, and knowing where to draw that line is the hard part.
- Avoid a standards process that slows implementation; with agents, specs
  and conformance tests might be a side effect of building ports.
- Treat the R implementation as the normative reference, and use the
  resulting specs and tests to assess other implementations.

The original posts are not reproduced here. The thread is on the Bioconductor
community Zulip; a ChatGPT conversation about shared cross-language test
fixtures also informed the E4 design.
