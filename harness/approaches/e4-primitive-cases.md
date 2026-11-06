E4 — SHARED CASES WITH PRIMITIVE INPUTS AND PROJECTED OUTPUTS.

- Conformance cases are data files (e.g. YAML) holding only primitive values:
  the inputs as constructor arguments (vectors of ints, strings, logicals),
  the operation and its arguments, and the expected result as a projection
  onto primitive vectors (e.g. start/end/seqnames/strand, hit pairs, mcols
  columns), plus a comparison mode per field (ordered, set, implementation-
  defined, ...). Objects themselves are never serialized.
- An R runner builds objects natively, runs the operation, projects the
  result, and fills in expected values (R is normative). A Python runner
  builds the Python objects, runs, projects, and checks.
- Each case carries a clause id; short ADRs record cases where a comparison
  mode other than exact is chosen or Python deliberately diverges.
- Implement the Python package until it passes the cases, alongside its own
  pytest suite.
