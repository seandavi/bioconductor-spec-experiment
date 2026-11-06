# Maintenance replay

Simulated upstream behavior change: `distance()` reverts to the pre-Bioc-2.13
rule (overlapping → 0; otherwise gap + 1, so adjacent → 1). A real historical
NEWS change, applied in reverse, so ports see "upstream changed".

- `Rprofile.replay`: patches the IntegerRanges `distance` method in `.First`;
  GRanges delegates to it. Versions are unchanged (a real release would bump them).
- `Rscript`: wrapper that sets `R_PROFILE_USER` (absolute paths for this machine).
- `detection/`: each approach's own drift check, run under the patched R.
