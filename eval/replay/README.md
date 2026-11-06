# Maintenance replay scoring

Upstream change simulated by `harness/replay/Rprofile.replay` (on main): `distance()`
reverts to the pre-Bioc-2.13 rule. Oracle re-run `--fresh` under that R: 93 expected
values change (all `distance`), 0 inputs change.

Ports scored here are the `replay/<exp>` branches (each port updated by its own
process from one release note); `ports.replay.toml` is the registry used. BiocPy is
unchanged (control). `results/` is the full scored output under the patched R.

| port | distance cases pass, before → after | regressions elsewhere |
|---|---|---|
| e0 | 89 → 89 / 93 | 0 |
| e0b | 71 → 71 / 93 | 0 |
| e1 | 93 → 93 / 93 | 0 |
| e2 | 88 → 88 / 93 | 0 |
| e3 | 84 → 84 / 93 | 0 |
| e4 | 84 → 84 / 93 | 0 |
| biocpy (not updated) | 37 → 0 / 93 | 0 |
