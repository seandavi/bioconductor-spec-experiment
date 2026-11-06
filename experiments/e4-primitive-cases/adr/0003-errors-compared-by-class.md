# ADR-0003: errors are compared by "raised", not by message

- Status: accepted.
- Golden `kind: error` means Python must raise `ValueError` or `TypeError`
  while building inputs or running the op. Messages are informational.
- R warnings are recorded in goldens (`warnings`) but never compared.
- Any other exception type (KeyError, AttributeError, ...) fails the case.
