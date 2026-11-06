E3 — SPEC AND ADRs FIRST, IMPLEMENTATION GUIDED BY THEM.

- First write a behavioral specification for the target classes and
  operations: language-independent, describing semantics (not R or Python
  code), with stable clause ids. Behavior the spec does not describe is, by
  definition, unimportant and implementation-defined.
- Record non-obvious semantic choices as ADRs (context, options, decision),
  each linked to clause ids.
- You may run R to check what the reference does while writing the spec.
- Then implement the Python package and its pytest suite, guided by the spec:
  tests cite the clause ids they exercise. The spec is not a test suite, and
  no cross-language harness or golden files are required.
