# Writing an adapter

An adapter translates eval calls (R vocabulary) into one port's public API and
projects the port's results back. It must be **thin**: if a result is wrong,
the port is wrong, and the score must say so.

## You get

| Read | Don't read |
|---|---|
| this file, `spec/PROJECTION.md`, `spec/ops.json`, `spec/CASES.md` | `build/expected.jsonl`, scored results, the report |
| the port's public API: README, docstrings, public names in `src/` | the port's tests, cases, goldens, ADR test pins |
| `build/resolved.jsonl` (inputs + calls, no expected values) | other adapters; `background/`; other experiments |
| output of `just adapter-check PORT` and `just roundtrip PORT` | |

You write one file: `eval/adapters/<port>.py`. Imports allowed: stdlib,
numpy, `adapter_base`, the port, and the port's own declared dependencies
(e.g. polars).

## Runner (already built; for orientation)

```
uv run --frozen --project <port dir> python eval/runner/run_port.py --port <port> ...
```

Per case: build each input fresh (`build`) → evaluate the call tree
(`$call` innermost first, via `call`) → `project` the result → one JSONL record.
The runner holds the 30 s timeout (`SIGALRM`, a `BaseException`, so a port's broad `except Exception` cannot swallow it).
An argument list R would reject (`BindError`) is a `port_error` with `origin: bind`; the oracle shows R errors there too.

## Interface

```python
from adapter_base import AdapterBase, RVector, Unsupported, UnknownArgument, MISSING

class Adapter(AdapterBase):
    port = "e1"                          # registry name
    port_packages = ("ranges_e1",)       # frames under these = "port raised"
    unsupported_ops = {"merge": "no public Seqinfo merge"}   # op -> reason, static

    # inputs: one per kind you support (vector -> RVector by default)
    def build_IRanges(self, p: dict): ...
    def build_GRanges(self, p: dict): ...
    def build_Seqinfo(self, p: dict): ...
    def build_DataFrame(self, p: dict): ...
    def build_Factor(self, p: dict): ...
    def build_IntegerList(self, p: dict): ...

    # ops: op_<ops.json method>, e.g. "[" -> op_subset, "names<-" -> op_set_names
    def op_findOverlaps(self, *, query, subject=MISSING, **kw): ...
    def op_c(self, *dots, x, **kw): ...

    def project(self, value) -> dict: ...  # native -> PROJECTION.md shape
```

## What an op method receives

The base binds the case's call against `ops.json` before calling you:

| Case has | You get |
|---|---|
| positional arg *k* | keyword `positional[k]` (R name, `.` → `_`: `ignore.strand` → `ignore_strand`) |
| kwarg in `named` | keyword, same renaming |
| extra positionals, op `dots` ≠ null | `*dots` (in order) |
| extra named args, `dots: "mcols"` | keyword `_mcols`: list of `(R name, value)` |
| anything else | never reaches you: R would reject it (`unused argument`); recorded as a port error |

Only arguments present in the case are passed. **Never fill in R defaults**:
the port's own defaults are under test. Values are:

| Value | Type |
|---|---|
| input / nested result | your native object (what `build_*` / an `op_*` returned) |
| atomic (incl. bare scalars) | `RVector(type, values, names)`; `.scalar()`, `.list()`, `.np()`, `.np_masked()` |
| R `NULL` | `None` |

## Projecting

`project` must handle every native the port can return. Use the base
helpers; they enforce PROJECTION.md (NA → null, `%.17g`, maximal Rle runs,
zero-length names → null):

| Helper | Use |
|---|---|
| `p_vector_from(arr_or_list, names=None)` | dtype → `type`: int* → int, float* → dbl, bool → lgl, str/object-of-str → chr; masked/None/NaN-in-int → null |
| `p_vector(type, values, names=None)` | explicit type |
| `p_null()` | |
| `p_iranges(start, width, names, mcols)` / `p_granges(seqnames, start, width, strand, names, mcols, seqinfo)` | `mcols` = list of `(name, column projection)` |
| `p_seqinfo(seqnames, seqlengths, is_circular, genome)` | |
| `p_dataframe(nrow, columns, rownames=None)` | `columns` = list of `(name, column projection)` |
| `p_factor(levels, codes, base=1)` | `base=0` for 0-based codes (negative = NA); output is 1-based |
| `p_integerlist(elements, names=None, shift=0)` | `shift` is added to every non-NA element (use `shift=1` for 0-based index lists such as `revmap`) |
| `p_hits(from_, to, n_left, n_right)` | **1-based** in |
| `p_rle(values, lengths=None, type=None)` / `p_rlelist(names, rles)` | dense or run form; normalized to maximal runs (NA and NaN each form their own runs) |

`RVector.np()` returns int64 / float64 / bool / object arrays and raises on NA in int/lgl/chr (dbl NA and NaN both become `nan`);
use `.np_masked()` when NA must survive. `RVector` and `None` results are projected by the base; `project` only sees port natives.

## Allowed vs forbidden

| Allowed (representation) | Forbidden (semantics) |
|---|---|
| R name / enum → port name / enum | branching on case id, source, tags or input *values* |
| index base shift where `ops.json` marks `index_valued`, per the port's documented base | shifting anything else (e.g. ranks) |
| `RVector` → list / ndarray / masked array; port containers → projection | coercing types toward R (float → int, str → factor) |
| dense coverage → runs; categorical → codes; mask/None → null | sorting, deduplicating, reordering, renaming port output |
| copy an input before an in-place port setter so the op returns a new object | computing a result the port did not compute |
| two public port calls when the port splits one R op (comment it) | catching a port exception to return a value; raising to imitate an R error; retrying with other args |
| `raise Unsupported(reason)` for an op or argument form the public API lacks | private/underscore port APIs (unless the port's README documents them) |

R subscripts (`[`): positive integers are shifted to the port's base.
Negative / zero / NA integer subscripts, logical recycling and name
subscripts are passed only if the port documents R semantics for them;
otherwise `Unsupported("subscript:<form>")`.

## Classification (runner, mechanical)

| Raised | Record status |
|---|---|
| `Unsupported(reason)`, or op listed in `unsupported_ops` | `unsupported` |
| `UnknownArgument(name)` from you (an R argument you can't map), or your signature rejects a legal R arg | `unsupported` (`argument:<name>`) |
| any exception whose traceback enters `port_packages` | `port_error` (phase `build` → always `port_error_unexpected`; see COMPARATOR.md) |
| any other exception (TypeError in translation, AttributeError on a missing port API, bad projection) | `adapter_error` |
| per-case timeout (default 30 s) | `port_error` (`origin: timeout`) |

`unsupported` and `adapter_error` never score as pass. Python warnings are
recorded, never compared.

## Dev loop and done

1. `just roundtrip PORT`: `project(build(p)) == p` for every resolved input.
2. `just adapter-check PORT`: lists only `adapter_error` and `unsupported` (no R values).
3. Done when: every `ops.json` op is an `op_*` or in `unsupported_ops`; round
   trip is clean for every supported input kind; zero `adapter_error` you can
   fix without semantic code; file header lists the port docs used, the
   port's index base, and every non-obvious translation.
4. Review gate: a reviewer audits the diff against the table above.

## Self-test adapters (eval-owned, not ports)

| Adapter | Behavior | Must |
|---|---|---|
| `r_echo` | overrides `run_case`: returns R's own projection from `expected.jsonl`; R error → `port_error` | score 100 % pass |
| `mutant` | `r_echo` + a deterministic perturbation on cases where `int(hash[:8], 16) % 5 == 0`; writes the predicted outcome per case | scored outcome == predicted, every case |

Mutant perturbations (first applicable in a hash-seeded cycle):

| Perturbation | Applies to | Predicted |
|---|---|---|
| bump first leaf (+1 / flip / append `"x"`) | R value with a leaf | `mismatch` |
| drop last element | R value, non-empty | `mismatch` |
| reverse `to` within one query | R2 Hits, a query with ≥ 2 hits | `pass` (R2 control) |
| swap to another universe hit | R3 case, a query with ≥ 2 universe hits | `pass` (R3 control) |
| return `NULL` | R error | `mismatch` |
| raise inside `fakeport` during call | R value | `port_error_unexpected` |
| raise inside `fakeport` during build | R error | `port_error_unexpected` |
| raise `TypeError` in adapter code | any | `adapter_error` |
| raise `Unsupported` | any | `unsupported` |

`fakeport` (`eval/selftest/fakeport/`) is a stdlib package that only raises,
so the traceback rule is exercised for real. Both adapters read `expected.jsonl` through `setup(ctx)` (passed only to
self-test ports) and override `run_case(case, stage)`; `r_echo` raises an R error through `fakeport`, so it is a `port_error`.
`mutant` writes `build/results/mutant.predicted.jsonl` (`id`, `predicted`, `perturbation`). Bump and drop are not applied to
R3 cases (any in-universe value passes, so `mismatch` is not guaranteed); `swap_universe` is their control.
`just selftest` fails unless `r_echo` scores 100 % and every mutant outcome equals its prediction.
