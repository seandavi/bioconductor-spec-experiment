E1 — PORT + ADVERSARIAL DIFFERENTIAL TESTING, R NORMATIVE.

- Port first; process is lightweight. Alongside the port, build a
  differential harness of your own design: it runs the same case in R and
  Python and compares results. Its form (case format, how values cross the
  language boundary, comparison rules) is yours to choose.
- Seed it with cases of your own. After implementation, a separate red-team
  agent will extend it with cases designed to break the port; design the
  harness so someone else can add cases easily.
- Divergences are resolved by fixing Python to match R, or, when matching R
  is wrong or not worth it, by a short ADR in adrs/ (context, options,
  decision). ADRs exist only where a real decision was made.
