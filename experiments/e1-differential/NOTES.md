# NOTES — e1-differential

Append-only. Timestamped judgment calls and human interventions.

- 2026-11-05 23:30 plan: harness cases in YAML (pyyaml dev dep); Python
  converts to JSON for R (jsonlite). R never parses YAML.
- 2026-11-05 23:30 plan: canonical compare form uses R's 1-based indices;
  adapter post-hooks convert Python results. R encoder stays generic.
- 2026-11-05 23:30 plan: R results cached in committed `r_cache/`, keyed by
  case hash + package versions, so `uv run pytest` runs without R and
  version bumps surface as cache diffs.
- 2026-11-05 23:30 plan: R probes — `reduce` keeps a zero-width range
  outside others (`20-19`) and absorbs one inside; zero-width `5-4` overlaps
  `1-10`, `11-10` does not; `gaps(GRanges)` with seqlengths emits full-length
  gaps on empty strands.
- 2026-11-05 23:55 harness: R-side warning handler records but must NOT muffle.
  Muffling changed R's control flow (`IRanges(start=2147483647, width=2)` with
  doubles errors only if the coercion warning propagates).
- 2026-11-05 23:55 R quirk: `IRanges(start=<int vec>)` alone builds runs of
  consecutive +1 integers (and TRUE-runs for logical); implemented, not an ADR.
- 2026-11-05 23:58 retracted: the apparent integer-overflow divergence was a
  harness artifact (warning muffling). R errors for integer inputs too; no ADR.
  Lesson: a divergence is only real after checking the harness did not alter R.
- 2026-11-06 00:20 R quirk (matched, no ADR): `reduce(min.gapwidth=0)` merges
  identical empty ranges but silently drops the repeats from `revmap`
  (e.g. 5/0,5/0,5/3 -> revmap [1],[3]). Modelled as "drop output range equal to previous".
  `gaps` and `disjoin` also ignore empty ranges as coverage; `disjoin` still
  breaks other ranges at an empty range's position (probed).
- 2026-11-06 01:10 overlaps: predicates fitted to R on an exhaustive 24x24 pool
  (types x maxgap x minoverlap), see `_kernels/overlaps.py` docstring. Facts found:
  `within` ignores maxgap<=0; minoverlap only applies when >0; `any` uses the
  two-sided gap test (so empty ranges behave); `any` with maxgap>=0 and minoverlap>0 errors.
- 2026-11-06 01:10 overlaps: hit order within a query is engine-dependent. R builds an
  NCList on the shorter side (ties -> subject). nq>=ns: subjects in (start asc, width desc,
  idx) order; nq<ns: index order. `select="arbitrary"` = first hit in that order
  (subject mode) or last-seen = max index (query mode). Matched exactly; no ADR needed.
  `first`/`last` are min/max subject index in both modes.
- 2026-11-06 01:10 ADR-0003: counts are plain arrays; R's query names stripped by the adapter.
- 2026-11-06 01:10 drop.redundant: modelled as "drop (i,j), i>j, if (j,i) is also a hit";
  matches R on symmetric (any/equal) and asymmetric (within/start) grids.
- 2026-11-06 01:50 nearest/precede/follow: ported from R source of the S4 methods (read via
  `selectMethod`): nearest uses findOverlaps(maxgap=0) first (adjacent counts), then
  precede(first)/follow(last); a tie goes to the follow side (arbitrary) or keeps both (all).
  Self mode (`nearest(x)`) uses findOverlaps(drop.self) whose select="arbitrary" is the max
  subject index (not the engine's first). distance() overflow -> NA (masked) + warning, as R.
- 2026-11-06 ordering: order/rank/duplicated use key (start, width), stable; `decreasing` keeps
  ties in index order (R radix). R's `rank(<IRanges>)` default errors -> ADR-0004. Not ported:
  multi-argument `order(x, y)`, `rank(ties.method="random")`, `==`/`<`/`match`/`xtfrm` (outside TARGET).
- 2026-11-05 23:53 correction: NOTES timestamps from the stage-4 entry to the ordering entry
  were estimated and ran ahead of the wall clock (real time at end of stage 9: 23:53, started
  23:26). Later entries use `date`.
- 2026-11-06 00:35 argument validation: a generated "bad args" grid (every op x NA/str/len2/neg/
  double/lgl/empty) showed 200/616 mismatches, almost all validation parity. Added `_args.py`
  (isTRUEorFALSE, isSingleNumber, match.arg with partial matching) + `@flags` decorator; also
  logical subscripts recycle, logical vectors are rejected where R wants numbers. 0 left except
  int32-overflow territory (ADR-0005 generalised; 2147483647 dropped from the generator).
- 2026-11-06 00:35 deviation from PLAN: psetdiff returns IRanges (R does too), so no list-return
  ADR; `compare: arbitrary` machinery exists in the harness but no case needs it (engine pick is modelled).
- 2026-11-06 00:50 GRanges parallel ops (found by a random seqinfo-pair grid, 640 cases): pintersect
  returns x's seqlevels only (values merged from y; y just validates/warns); pintersect strand =
  merged only for hits, else x's; psetdiff strand = merged always; punion strand = x's.
- 2026-11-06T05:46Z red team: baseline 3808 passed. Added 33 circular-seqinfo cases;
  17 expose lost periodic/origin hits (including nearest and empty ranges), 16 match R.
  Direct R execution confirms origin overlap. Automatic circular-seqinfo semantics are
  in scope; explicit circle.length arguments remain excluded. These are bugs, not ADRs.
- 2026-11-06T05:46Z red team: unresolved cases use reasoned, strict xfail marks via
  a collection hook. They still execute both sides; fixes produce XPASS failures.
  Existing expectations and src remain unchanged. No divergences were accepted.
- 2026-11-06T05:44:49Z correction: the two preceding red-team entries used an estimated
  05:46 timestamp; the first red-team commit actually completed at 05:42:33Z.
- 2026-11-06T05:44:49Z red team: 40 metadata/validation cases: 25 strict xfails, 15 passes.
  Found metadata recycling errors, overwritten-hit collision, lost masked logical NA,
  coverage column/NA handling gaps, invalid Seqinfo/strand coercions, and unrecycled names.
  Seqinfo replacement silently remaps/drops levels: judgment call requiring an ADR,
  not an accepted divergence. Numeric-name case uses direct constructors to avoid
  the encoder's unrelated inability to JSON-serialize numpy scalar sequence names.
- 2026-11-06T05:48:41Z red team: 46 empty/zero-width/ordering/negative-coordinate cases;
  36 pass, 10 strict xfail. Default constructors, blank-name indexing, and typed NA
  coverage weights diverge. Negative genomic nearest-family behavior needs an ADR:
  R loses hits or errors; Python gives geometric answers. No decision accepted.
- 2026-11-06T05:48:41Z harness correction: typed empty numpy metadata became Polars Null
  through tolist(), falsely appearing as character metadata. Preserve array dtype
  in the literal decoder; seven regression tests cover it. Four apparent port
  failures now pass. Xfail only catches comparison failures, not harness crashes.
- 2026-11-06T05:48:41Z validation: full suite 3888 passed, 52 xfailed. Source unchanged;
  existing test expectations unchanged. All 119 red-team cases ran in R and Python.
- 2026-11-06T05:52:43Z red-team report: per-case tables cover all 119 cases. Classified 47
  bug cases, 5 unresolved ADR judgments, 67 agreements. Renamed three trim cases
  to say circularity (their seqlength is known); inputs/results are unchanged.
  Fresh-R new-case run: 67 passed, 52 xfailed. Cache-only full run: 3888 passed,
  52 xfailed. Unmasked origin-overlap case fails as expected. Regeneration and
  compileall pass. No source edits or accepted ADRs; independent review pending.

## Review round (independent review + red team)

- 2026-11-06 findings 1-12 handled; commits reference finding numbers. Tests: 4746 pass, 0 xfail.
- Finding 1 (cache fingerprint): one fingerprint = R versions + digest of the R scripts, persisted and compared;
  offline mode checks the script digest. Unit tests in `tests/unit/test_cache_lifecycle.py`.
- Finding 2 (immutability): arrays are backed by immutable `bytes` buffers; accessors return unwritable views; mcols are cloned at
  every public boundary. Tests in `tests/unit/test_immutability.py`.
- Findings 4/5/8/9/10/11/7: short mcols recycle (warning when not a multiple), recycled GRanges coordinates recycle names, flank
  rejects masked logicals, literal strand "NA"/numeric circularity/numeric Seqinfo names are errors, `IRanges()`/`GRanges()` are
  empty, blank names never match a subscript, pintersect replaces an existing `hit` column.
- Finding 6 (coverage NA): reviewer suggested local NA effects; R does not do that. Reverse-engineered on 4000 random layouts:
  R has two engines ("sort" iff width >= 4 * length(x), else "hash") and NA stays local only if the ranges form a partition of
  1..width, otherwise everything from the first NA range's start is NA. Implemented exactly (0/4000 mismatches) rather than the
  reviewer's local-NA suggestion; it is what R does.
- Finding 3 (circular): hits on any copy s + kL after moving both sides to the first circle; circular hits are ordered by
  subject index; `equal` is start and end modulo L independently; `arbitrary` = first found by (copy order 0,-1,+1,.., then
  NCList order), or last-seen when the query side is shorter. Remaining engine-path-dependent maxgap behaviour: ADR-0008.
- Finding 12: `seqinfo<-` now follows R (levels may only be appended; no ADR needed). Negative-coordinate nearest family: ADR-0007.
- Harness bugs found while fixing: (a) R literal decoder always passed strand "*" so every empty GRanges input errored in R and
  ~30 cases were vacuous (both sides erred); (b) `grid_seqinfo` called nearest with select="first" (R error, 30 vacuous cases).
  Both fixed; GRanges precede/follow with an empty side return empty Hits and nearest an all-NA logical vector, as R.
- Rejected: "four historical stage commits omit the attribution trailer" -- history may not be rewritten (PROTOCOL: no amend, no
  rebase); the omission is in commits already shared. Trailers are present on every commit of this round.
- Rejected: reviewer's "do not blindly reproduce sentinel artifacts" is accepted for cross-sequence hits and errors (ADR-0007), but
  the narrow "coordinate 0 on a non-first level drops the hit" quirk stays because ordinary data (empty range at position 1) triggers it.

## Maintenance replay: IRanges 2.47.3 / GenomicRanges 1.65.2 (distance)

- 2026-11-06T16:22Z upstream change: `distance()` is 0 for overlapping ranges, otherwise gap width + 1 (adjacent = 1);
  GRanges delegates to the IntegerRanges method; nearest/precede/follow unchanged. Followed README "Maintaining": reran R,
  read `git diff r_cache/`, fixed Python. `--refresh-r` changed 40 cases, all `distance` (IRanges + GRanges, incl. NA/strand
  cases and the nearest/precede grids that call it); nothing else moved. Python: `d = max(max(start) - min(end), 0)`. No ADR:
  Python matches R. Full suite 4752 passed (4746 + 6 new).
- Harness quirk: the R packages on this machine already carry the new behaviour but still report IRanges 2.47.2 / GenomicRanges
  1.65.1, so the cache fingerprint did not change and plain `uv run pytest` stayed green on stale R results. Only
  `--refresh-r` exposed the change. A version key alone cannot catch behaviour changes shipped under an unchanged version.
- Case id `distance-adjacent-is-zero` (nearest.yaml) now expects 1; the id is kept so the cache key/history stay stable (cases are
  append-only). The new `distance-adjacent-is-one` and friends cover the new rule.
- Found: `distance-basic` and `distance-zero-width` had `width:` outside the `IRanges` literal (R ignored it); fixed. The same
  slip remains in nearest-far-both-sides, nearest-zero-width-subject, nearest-unsorted-subject, precede-zero-width-x
  (nearest.yaml) and disjoin-zero-width-only (inter.yaml): not touched here, they are unrelated to distance and still need review.
- PLAN.md line 106 ("distance = gap length, 0 if overlapping/adjacent") is the historical plan; left as written.
