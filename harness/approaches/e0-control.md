E0 — CONTROL: a plain agentic port, the way a competent developer would do it
today without any cross-language process.

- Build the Python package and its own pytest suite.
- You may consult R documentation (`Rscript -e '?reduce'`, help pages, package
  source under the R library) and run small ad hoc R snippets to clarify what a
  function does, as any porter would.
- Do NOT build any systematic R-vs-Python comparison: no harness that runs both,
  no shared case files, no golden files generated from R, no formal spec or ADR
  process. Expected values in tests are written by you.
- Ordinary design notes in README/NOTES are fine.
