"""Adjudication verdicts: one rule per cluster signature (regex, first match wins).

Each rule: (signature regex, class, evidence, rationale). Classes:
port_bug | declared_divergence | r_accident | adapter_fault | harness_fault.
Doc citations are `<port>:<path>:<line>` relative to the port's project dir (eval/ports.toml).
R evidence names an R probe in this directory or an inline `Rscript` result (R 4.6 / IRanges 2.47.2 /
GenomicRanges 1.65.1).
"""
from __future__ import annotations

PB, DD, RA, AF, HF = "port_bug", "declared_divergence", "r_accident", "adapter_fault", "harness_fault"

# R facts reused across ports (each checked in R for this adjudication)
R_COV_NA = ("probe_coverage_na.py: in 228/372 NA-weight cases R puts NA where no NA-weighted range covers "
            "(IRanges(c(7,7,6),width=c(1,1,5)), weight=c(2,NA,3) -> pos 8-10 NA, only covered by weight 3); "
            "swapping the first two ranges of a partition changes [NA,3,1,NA,4] to all-NA")
R_RANK = ("rank(IRanges(c(5,1,5),width=c(2,3,2))) errors \"'arg' must be of length 1\": S4Vectors rank,Vector "
          "passes the whole ties.method choice vector to .local; base rank() defaults to 'average'")
R_SENTINEL = ("GRanges precede/follow/nearest put all seqlevels on one line with sentinels; coordinates <= 0 "
              "cross them: heldout/nearest/GRanges/76 returns a chrX query's nearest on chrY, "
              "heldout/nearest/GRanges/110 a chrY query's nearest on chr3")
R_INT32 = ("R does coordinate arithmetic in int32: disjoin errors \"'start(q)' and 'end(q)' must have the same "
           "length\" on valid input ending at 2^31-1; gaps/setdiff lose the hole; findOverlaps(maxgap=2^31-1) "
           "finds nothing (e4:adr/0010-r-overflow-defects-not-ported.md)")
R_RESIZE_FIX = ("resize(GRanges('chr1', IRanges(c(28,14), width=c(8,12)), '+'), 3, fix=c('end','center'), "
                "ignore.strand=TRUE) starts 33 23; ignore.strand=FALSE (same, all '+') and IRanges resize give 33 18")
R_RESIZE_EMPTY = ("resize(<empty GRanges>, 12, fix=c('start','bogus')) returns; the same fix on a non-empty "
                  "GRanges errors: validation depends on length(x)")
R_INT32_NA = ("R wraps int32 overflow to NA and returns objects its own validity check rejects: "
              "validObject(promoters(IRanges(2147483642, width=3), 0, 100)) and flank(IRanges(2147483638, width=10), 5, "
              "start=FALSE) fail \"'start(x)', 'end(x)', and 'width(x)' cannot contain NAs\"; shift(GRanges('chr1:10-14'), "
              "2147483647) silently has start NA while the same IRanges shift errors")
R_EMPTY_TYPE = ("GRanges precede/follow return SortedByQueryHits when the query or subject is empty (integer "
                "otherwise; IRanges: integer either way) and GRanges nearest returns logical NA for an empty "
                "subject (IRanges: integer NA): R's return type depends on emptiness")
R_RESOLVE_EMPTY = ("pintersect(x, x, resolve.empty=NA|TRUE|'a'|1:2|integer(0)) returns; the same values error once an "
                   "empty intersection must be resolved (pintersect(x, y) with a non-overlapping pair)")
R_CIRCULAR_PROBE = ("probe_circular.py: with every isCircular set FALSE R returns exactly the port's value, i.e. the "
                    "port computes linear semantics on a circular sequence")
R_TRIM_NAMES = "trim(x, use.names=FALSE) keeps names(x) (R ignores the argument); shift(use.names=FALSE) drops them"

RULES: list[tuple[str, str, str, str]] = [
    # ------------------------------------------------------------------ cross-port R accidents
    (r"\|coverage\|.*\|r_nonlocal_na$", RA, R_COV_NA + "; the port's value is the local-NA sum with R's type",
     "R's NA-weight coverage is not a sum here (NA leaks past the NA range / depends on order)"),
    # R's value is defective AND the port has its own defect on the same case: the port's defect counts
    (r"^e2\|coverage\|.*\|r_nonlocal_na\+port_wrong$", DD,
     "e2:NOTES.md:27 'coverage: ... NA weight gives garbage (not pinned)': e2 refuses NA weights (IRanges) or computes "
     "with NaN (GRanges); " + R_COV_NA, "R's NA coverage is defective; e2's NA-weight handling is its declared one"),
    (r"\|coverage\|.*\|r_nonlocal_na\+port_wrong$", PB,
     "probe_coverage_na.local_na applied to the port's value: not the local-NA sum with R's type and lengths (e1: NaN "
     "in a double Rle where R has NA_integer_/NA_real_, same positions as R; e4/BiocPy: NA weights as 0 or a crash). "
     + R_COV_NA, "R's NA coverage is defective, but the port is wrong there too (crash / NaN / NA as 0 / type)"),
    (r"\|nearest\|.*\|r_translation_dependent\+port_wrong$", PB,
     "probe_translate.py r_shifted (+ select='all' universe for R3): the port's value differs from R's answer for the "
     "same case shifted by +1000, e.g. heldout/precede/GRanges/40 R 6 5 3 NA 2, shifted 6 5 4 8 2, e0b/BiocPy 2 5 3 8 9",
     "R's value is a sentinel artefact, but the port also misses the translation-invariant answer"),
    (r"\|nearest\|r_error_port_value\|na_in_logic\|coord_le0\+port_wrong$", PB,
     "e1/redteam_edges/rt-edge-gr-negative-follow: R errors; shifted +1000 R gives 2->1, 4->3 (as e0-e4); the port "
     "adds 3->2", "R fails on valid input, and the port misses R's shifted answer"),
    (r"^(e2|biocpy)\|overlaps\|mismatch\|diff\|maxgap_int32_max$", RA,
     "findOverlaps(IRanges(1,1), IRanges(2147483642,2147483642), maxgap=2147483645) hits; maxgap 2147483646 or "
     "2147483647 loses the hit (non-monotonic in maxgap: int32 overflow). The port finds the hit. " + R_INT32,
     "R's maxgap near 2^31 overflows"),
    (r"\|setops\|port_error:call\|rejects_valid_resolve_empty", PB,
     "R's resolve.empty accepts 'start.x' and its unique prefix 'start' (match.arg); e4 spells the choice 'start' and "
     "rejects 'start.x' (src/ranges_e4/setops.py:66-67, clauses.md:88), e0/e2/e3 reject the partial 'start'",
     "valid resolve.empty value rejected"),
    (r"\|nearest\|mismatch\|diff\|empty$", RA, R_EMPTY_TYPE,
     "empty query/subject: R's GRanges result class (Hits / logical) differs from its own non-empty and IRanges answers"),
    (r"^(e0|e0b|e2|e3|e4|biocpy)\|intra\|mismatch\|diff\|strand$", RA, R_RESIZE_FIX,
     "R's ignore.strand=TRUE resize ignores the per-range fix vector"),
    (r"^(e0|e0b|e2|e3|e4|biocpy)\|intra\|port_error:call\|ValueError:'?fix'? must", RA, R_RESIZE_EMPTY,
     "R skips fix validation on empty input only"),
    # ------------------------------------------------------------------ e1
    (r"^e1\|coverage\|mismatch\|nan_as_na\|na_weight_dbl$", PB,
     "e1:src/ranges_e1/iranges.py:512 weights become float NaN; Rle values keep NaN; README.md:30 says NA is a "
     "masked entry; e1's own harness hides it (tests/differential/adapter.py:241 'NaN = NA')",
     "NA coverage is returned as NaN (R: NA_real_) where R's value is the principled local-NA sum"),
    (r"^e1\|coverage\|mismatch\|(nan_as_na\+)?num_type\|(na_weight_int|-)$", PB,
     "e1:src/ranges_e1/iranges.py:526 converts integer weights to float64 whenever any weight is NA",
     "integer weights with any NA give a double Rle (R keeps integer, NA_integer_)"),
    (r"^e1\|coverage\|port_error:call\|ValueError:'width' must be NULL", PB,
     "R accepts coverage(width=1.5) (truncates) and positional width=NA (as NULL); no e1 doc",
     "rejects widths R accepts"),
    (r"^e1\|nearest\|mismatch\|diff\|(r_translation_dependent|coord_le0|int32_edge)$", DD,
     "e1:adrs/0007-nearest-family-nonpositive-coordinates.md; README.md:73. " + R_SENTINEL,
     "geometric nearest/precede/follow for coordinates <= 0 by ADR-0007"),
    (r"^e1\|nearest\|r_error_port_value\|na_in_logic", DD,
     "e1:adrs/0007-nearest-family-nonpositive-coordinates.md ('follow(select=\"all\") fails with \"logical "
     "subscript contains NAs\"')", "ADR-0007 returns the geometric answer where R errors"),
    (r"^e1\|overlaps\|mismatch\|names\|-$", DD, "e1:adrs/0003-plain-arrays-no-names.md; README.md:69",
     "countOverlaps returns unnamed counts by ADR-0003"),
    (r"^e1\|overlaps\|mismatch\|diff\|circular$", DD,
     "e1:adrs/0008-circular-maxgap-hits.md (pins div-circ-end-maxgap-*); R: same pair hits alone, misses in context",
     "circular maxgap>0 end/start: one pairwise rule where R is engine-path dependent"),
    (r"^e1\|(overlaps|setops|inter)\|mismatch\|diff\|(int32_edge|maxgap_int32_max)$", DD, "e1:adrs/0005-coordinate-overflow-is-an-error.md; " + R_INT32,
     "int64 arithmetic near 2^31 by ADR-0005 (R wraps)"),
    (r"^e1\|\*\|r_error_port_value\|length_validation", DD, "e1:adrs/0005-coordinate-overflow-is-an-error.md; " + R_INT32,
     "disjoin at end=2^31-1: R's int32 end+1 overflow errors; e1 computes exactly (ADR-0005)"),
    (r"^e1\|intra\|port_error:call\|int32_overflow", DD, "e1:adrs/0005-coordinate-overflow-is-an-error.md; README.md:71",
     "raises where R wraps coordinates to NA (ADR-0005)"),
    (r"^e1\|setops\|mismatch\|diff\|out_of_bounds$", PB,
     "R GRanges intersect/setdiff = gaps(union(gaps(x),y)) bounded by [1, seqlength] (e2:NOTES.md GR.SETOPS.CLIP); "
     "e1 has no clipping in its set-op path and no ADR",
     "GRanges intersect/setdiff do not clip to [1, seqlength] as R does"),
    (r"^e1\|setops\|mismatch\|names\|-$", PB, "R: p* ops on unnamed x take y's names; no e1 doc",
     "pintersect/punion/psetdiff drop y's names when x has none"),
    (r"^e1\|ordering\|r_error_port_value\|rank_default_ties", DD, "e1:adrs/0004-rank-default-ties-method.md; " + R_RANK,
     "rank() defaults to 'average' by ADR-0004"),
    (r"^e1\|accessors\|mismatch\|names\|-$", DD, "e1:README.md:86 ('Not implemented: ... mcols row names'); "
     "adrs/0006-mcols-with-na-names.md:15", "mcols(x) has no rownames (R: names(x))"),
    (r"^e1\|accessors\|r_error_port_value\|na_input_rejected", DD, "e1:adrs/0006-mcols-with-na-names.md",
     "mcols() with NA names works by ADR-0006"),
    (r"^e1\|accessors\|mismatch\|diff\|(zero_width|circular)$", PB,
     "R: mcols(<GRanges without columns>) is a 0-column DataFrame with nrow = length(x) (IRanges: NULL); no e1 doc",
     "GRanges mcols() with no columns returns None (R: DataFrame)"),
    (r"^e1\|construct\|mismatch\|mcols\|-$", PB, "R: GRanges(seqnames, ranges=<IRanges with mcols>) keeps mcols(ranges)",
     "GRanges constructor drops the mcols of its `ranges` argument"),
    (r"^e1\|construct\|mismatch\|(seqinfo|seqlevel_order)\|-$", PB,
     "R: seqlevels follow names(seqlengths) or factor levels when given (e2/granges-core/gr.ctor.seqlengths, "
     "gr.ctor.factor-seqnames-levels); e1 uses first appearance in seqnames",
     "ignores the seqlevel order implied by seqlengths names / factor levels"),
    (r"^e1\|construct\|port_error:call\|AttributeError", PB, "R: GRanges(seqinfo=<character>) is accepted",
     "crashes on a character seqinfo"),
    (r"^e1\|construct\|port_error:call\|ValueError:the length of the supplied 'seqlengths'", PB,
     "R: GRanges('chr2', ..., seqlengths=c(chr2=151, chr3=197, chrY=152)) adds the extra levels",
     "rejects seqlengths naming more sequences than seqnames use"),
    (r"^e1\|construct\|r_error_port_value\|other:'seqlengths' must have names", PB,
     "R errors on unnamed seqlengths", "accepts unnamed seqlengths"),
    (r"^e1\|\*\|r_error_port_value\|type_validation", PB, "R: x[character(0)] errors when names(x) is NULL",
     "character subscript on an unnamed object returns empty instead of erroring"),
    (r"^e1\|seqinfo\|port_error:call\|ValueError:'genome' must be", PB, "R: genome(x) <- 1L stores \"1\"",
     "rejects a non-character genome R coerces"),
    (r"^e1\|intra\|mismatch\|diff\|strand$", RA, R_RESIZE_FIX,
     "R's ignore.strand=TRUE resize ignores the per-range fix vector; e1 matches the IRanges/strand-FALSE answer"),
    (r"^e1\|intra\|port_error:call\|ValueError:'fix' must be", RA, R_RESIZE_EMPTY,
     "R skips fix validation on empty input only"),
    (r"^e1\|seqinfo\|mismatch\|names\|-$", RA, R_TRIM_NAMES, "R ignores trim(use.names=FALSE); e1 honours it"),
    (r"^e1\|subset_concat\|mismatch\|mcol_types", PB,
     "R: c() of GRanges where only an empty object has column `flag` (chr(0)) gives a logical NA column",
     "fills a column from an empty object with its type (chr) instead of R's logical NA"),
    (r"^e1\|subset_concat\|port_error:call\|TypeError:concat needs", PB, "R: c(<IRanges>, <IRanges>, 1L) coerces 1L",
     "c() rejects an integer operand R coerces to IRanges"),

    # ------------------------------------------------------------------ e0
    (r"^e0\|(overlaps|nearest|coverage)\|mismatch\|equals_R_linear\|circular$", DD,
     "e0:README.md 'Not implemented: ... circular-sequence wrap-around'; NOTES.md 01:30 'circle.length/circular wrap'. "
     + R_CIRCULAR_PROBE, "no circular wrap-around (declared)"),
    (r"^e0\|(overlaps|nearest)\|mismatch\|diff\|circular$", DD,
     "e0:README.md 'Not implemented: ... circular-sequence wrap-around'; cases are circular overlaps whose "
     "select='arbitrary' pick or subset differs from both R and R-linear only by the circular hits",
     "no circular wrap-around (declared)"),
    (r"^e0\|\*\|r_error_port_value\|length_validation", DD,
     "e0:README.md 'Not implemented: circular-sequence wrap-around' (R errors on coverage width > circular length) and "
     "README.md Storage 'int64 ... R's overflow errors are not reproduced' (disjoin at 2^31-1)",
     "R errors from circular-length / int32 checks the port declares it does not implement"),
    (r"^e0\|coverage\|mismatch\|nan_as_na\+num_type\|na_weight_dbl$", PB,
     "e0:src/ranges_e0/_coverage.py:28 + GRanges per-seqlevel default: seqlevels without ranges get an int Rle "
     "although weights are double (R: double for every seqlevel)", "empty seqlevels of a double-weighted coverage are integer"),
    (r"^e0\|coverage\|mismatch\|num_type\|(-|na_weight_dbl)$", PB,
     "e0:src/ranges_e0/_coverage.py:28 ('Integral values give int64'): integral double weights and empty "
     "seqlevels give integer Rles (R keeps double)", "double weights produce integer coverage"),
    (r"^e0\|coverage\|mismatch\|diff\|na_arg$", PB,
     "R with NaN weights (the port's actual input) keeps NaN local: coverage(IRanges(c(1,2,5,8,11),width=c(1,3,3,3,2)), "
     "weight=c(NaN,3,1,NaN,4)) = NaN 3 1 NaN 4; e0 returns all-NaN (NaN enters the running sum and never leaves)",
     "a NaN weight poisons every later position"),
    (r"^e0\|coverage\|r_error_port_value\|other:", PB,
     "R rejects width=NA/TRUE, shift longer than x, and named width not identical to seqlevels(x); e0 accepts them",
     "missing argument validation in coverage"),
    (r"^e0\|\*\|r_error_port_value\|(flag_validation|scalar_validation|type_validation)", PB,
     "R validates flags (isTRUEorFALSE), maxgap/minoverlap/min.gapwidth scalars, weight/start/isCircular types; "
     "e0 coerces NA/strings/ints/negatives silently (no doc)", "missing argument validation"),
    (r"^e0\|(nearest|overlaps)\|r_error_port_value\|na_in_logic\|na_arg$", PB,
     "R errors 'missing value where TRUE/FALSE needed' on ignore.strand=NA / invert=NA; e0 accepts NA",
     "NA flag accepted"),
    (r"^e0\|\*\|port_error:call\|rejects_recycling", DD,
     "e0:README.md 'Only length-1 arguments are recycled (R also recycles other lengths with a warning)'; NOTES.md 23:20",
     "no cyclic recycling (declared)"),
    (r"^e0\|\*\|port_error:call\|rejects_names_padding", PB,
     "R pads short names with NA (IRanges(..., names='a') on 2 ranges, names(x) <- shorter); this is not recycling, "
     "so README's length-1 rule does not cover it", "names shorter than x are rejected instead of NA-padded"),
    (r"^e0\|\*\|port_error:call\|rejects_double_arg", DD,
     "e0:NOTES.md 23:20 'Fractional floats rejected'; R truncates (IRanges(start=1.5) -> 1) and accepts logical start runs",
     "fractional / logical coordinates rejected (declared)"),
    (r"^e0\|inter\|port_error:call\|TypeError:only N-dimensional", PB,
     "R recycles a gaps start/end vector over seqlevels; e0 crashes (TypeError) on a length-2 bound",
     "crash on per-seqlevel bound vectors"),
    (r"^e0\|nearest\|mismatch\|diff\|(r_translation_dependent|int32_edge|coord_le0)$", RA,
     "probe_translate.py: R's answer changes when every range is shifted +1000 (or coordinates near 2^31 overflow). "
     + R_SENTINEL, "R's GRanges nearest-family sentinel collisions; e0 gives the geometric answer"),
    (r"^e0\|nearest\|r_error_port_value\|na_in_logic\|coord_le0$", RA, R_SENTINEL + "; follow(select='all') errors "
     "'logical subscript contains NAs' on valid input", "R fails on valid input with coordinates <= 0"),
    (r"^e0\|overlaps\|mismatch\|names\|-$", PB, "R names countOverlaps by names(query); no e0 doc",
     "countOverlaps drops query names"),
    (r"^e0\|(overlaps|nearest|seqinfo)\|r_error_port_value\|seqinfo_incompatible", PB,
     "R errors on incompatible seqlengths/circularity and on seqinfo<- that drops/reorders levels; e0 NOTES only "
     "describes merge warnings", "merges incompatible Seqinfo / accepts invalid seqinfo<-"),
    (r"^e0\|seqinfo\|r_error_port_value\|other:The seqlevels", PB, "R refuses a reordered Seqinfo in seqinfo<-",
     "accepts a reordered Seqinfo"),
    (r"^e0\|setops\|mismatch\|seqinfo\|-$", PB,
     "R: GRanges pintersect keeps x's seqlevels only (e1:NOTES.md 00:50); e0 merges y's levels into the result",
     "pintersect result seqinfo includes y's extra seqlevels"),
    (r"^e0\|setops\|mismatch\|diff\|out_of_bounds$", PB,
     "R GRanges intersect/setdiff clip to [1, seqlength] (gaps-based definition); e0 does not clip", "no clipping in set ops"),
    (r"^e0\|setops\|mismatch\|names\|-$", PB, "R: punion on unnamed x takes y's names", "p* op drops y's names"),
    (r"^(e0|e0b)\|setops\|port_error:call\|(resolve_empty_validation|ValueError:The truth value)", RA, R_RESOLVE_EMPTY,
     "R validates resolve.empty only when an empty intersection occurs"),
    (r"^e0\|(setops|inter|overlaps|nearest|intra)\|mismatch\|diff\|(int32_edge(_distance)?|maxgap_int32_max)$", DD,
     "e0:README.md Storage 'int64 (R: int32, so R's overflow errors are not reproduced)'. " + R_INT32,
     "int64 arithmetic near 2^31 (declared)"),
    (r"^e0\|(intra|construct)\|r_error_port_value\|int32_overflow", DD,
     "e0:README.md Storage 'int64 ... R's overflow errors are not reproduced'", "no int32 overflow error (declared)"),
    (r"^e0\|subset_concat\|mismatch\|mcols\|-$", PB,
     "R coerces a double column to character with 15 significant digits (as.character(1234567.123456789) = "
     "'1234567.12345679'); e0 uses Python repr", "c() of chr+dbl columns formats doubles differently"),
    (r"^e0\|subset_concat\|mismatch\|mcol_types", PB,
     "R: c() where only an empty object has column `flag` gives logical NA", "takes the empty object's column type"),
    (r"^e0\|subset_concat\|port_error:call\|TypeError", PB, "R: c(<IRanges>, <IRanges>, 1L) coerces 1L",
     "c() rejects an integer operand"),
    (r"^e0\|metadata\|mismatch\|names\|-$", PB, "R keeps NA names; e0 stores NA names as '' (adapters/e0.md notes it, "
     "port docs do not)", "NA names become empty strings"),
    (r"^e0\|construct\|mismatch\|names\|-$", PB,
     "R: NA names kept, logical names become 'TRUE'/'FALSE', a recycled length-1 ranges keeps its name; e0: '', "
     "'True', dropped", "constructor name handling"),
    (r"^e0\|accessors\|mismatch\|names\|-$", PB, "R mcols(x) has rownames = names(x); no e0 doc",
     "mcols() without rownames"),
    (r"^e0\|accessors\|mismatch\|diff\|(zero_width|circular)$", PB,
     "R: mcols(<GRanges without columns>) is a 0-column DataFrame with nrow = length(x)", "returns NULL"),
    (r"^e0\|construct\|mismatch\|mcols\|-$", PB, "R: GRanges(seqnames, ranges=<IRanges with mcols>) keeps mcols(ranges)",
     "drops the mcols of `ranges`"),
    (r"^e0\|construct\|mismatch\|seqlevel_order\|-$", PB, "R uses factor levels as seqlevels (gr.ctor.factor-seqnames-levels)",
     "ignores factor level order"),
    (r"^e0\|construct\|mismatch\|diff\|plain$", PB,
     "R: IRanges(start=<int>) alone is the integer-set form (runs: c(1,2,3,7,8) -> 1-3, 7-8); e0:NOTES.md 23:20 "
     "says 'IRanges(5) (start only) = width 1 (R behaviour)', true only for length 1", "start-only constructor not run-encoded"),
    (r"^e0\|construct\|port_error:call\|ValueError:at least two", PB, "R: IRanges() is the empty IRanges",
     "IRanges() with no arguments raises"),
    (r"^e0\|construct\|port_error:call\|ValueError:strand values", PB, "R: NA strand becomes '*'", "NA strand rejected"),
    (r"^e0\|construct\|port_error:call\|AttributeError", PB, "R: GRanges(seqinfo=<character>) is accepted",
     "crashes on a character seqinfo"),
    (r"^e0\|construct\|r_error_port_value\|na_input_rejected", PB, "R errors on NA seqnames / NA Seqinfo names",
     "accepts NA seqnames"),
    (r"^e0\|ordering\|mismatch\|num_type\|-$", DD,
     "e0:README.md 'rank defaults to ties_method=\"first\"' (integer); R's GRanges rank default is 'average' (double)",
     "rank default ties method (declared)"),
    (r"^e0\|ordering\|r_error_port_value\|rank_default_ties", DD, "e0:README.md 'rank defaults to ties_method=\"first\"'; "
     + R_RANK, "rank() works without ties.method (declared)"),
    (r"^e0\|(nearest|ordering|overlaps)\|port_error:call\|ValueError:'(select|ties_method|type)' must be one of", PB,
     "R match.arg partial-matches 'a' (select -> 'all', type -> 'any', ties.method -> 'average')",
     "no partial matching of choice arguments"),
    (r"^e0\|overlaps\|mismatch\|diff\|plain$", PB, "R truncates a fractional minoverlap (5.5 -> 5); e0 compares with 5.5",
     "fractional minoverlap not truncated"),
    (r"^e0\|overlaps\|mismatch\|diff\|(coord_le0_r_invariant|zero_width)$", PB,
     "probe_translate.py: R's answer is translation-invariant; e0 differs on findOverlaps(type='within', drop.self, "
     "drop.redundant) with zero-width ranges", "drop.redundant/drop.self with zero-width ranges"),
    (r"^e0\|inter\|mismatch\|diff\|drop_empty_ranges$", PB,
     "R: with drop.empty.ranges=TRUE empty ranges still bridge/extend clusters before being dropped "
     "(e1:README.md 'R quirks reproduced'; e3:adr/0007); e0 drops them first", "empty ranges do not bridge in reduce(drop.empty.ranges=TRUE)"),
    (r"^e0\|inter\|mismatch\|diff\|min_gapwidth_0$", PB,
     "R reduce(min.gapwidth=0) merges identical empty ranges (e1:NOTES.md 00:20, e3:adr/0013); e0 keeps duplicates",
     "duplicate empty ranges not merged at min.gapwidth=0"),
    (r"^e0\|intra\|r_error_port_value\|other:'length\(shift\)'", PB, "R errors when shift is longer than x (also for empty x)",
     "no length check"),

    # ------------------------------------------------------------------ e0b
    (r"^e0b\|(overlaps|nearest|coverage)\|mismatch\|equals_R_linear\|circular$", DD,
     "e0b:README.md:84,102 'Circular wraparound is not implemented'; NOTES.md:27. " + R_CIRCULAR_PROBE,
     "no circular wrap-around (declared)"),
    (r"^e0b\|(overlaps|nearest|coverage|inter)\|mismatch\|diff\|circular$", DD,
     "e0b:README.md:84,102 'Circular wraparound is not implemented'; examined cases differ by circular hits / "
     "arbitrary picks among them, or gaps on circular sequences", "no circular wrap-around (declared)"),
    (r"^e0b\|\*\|r_error_port_value\|length_validation\|circular$", DD,
     "e0b:README.md:102 (no circular wrap-around): R's 'width cannot be greater than length of underlying circular "
     "sequence' check does not exist in the port", "circular check not implemented (declared)"),
    (r"^e0b\|\*\|r_error_port_value\|length_validation\|int32_edge$", DD, "e0b:README.md:52 'int64, no overflow check'. "
     + R_INT32, "int64 arithmetic (declared)"),
    (r"^e0b\|\*\|r_error_port_value\|length_validation\|-$", PB,
     "R errors when Seqinfo seqlengths/isCircular are shorter than seqnames; e0b recycles them", "lenient Seqinfo lengths"),
    (r"^e0b\|\*\|port_error:call\|rejects_recycling", DD,
     "e0b:README.md:103 'Scalar-or-vector arguments are broadcast with length 1 or n only (no R-style full recycling)'",
     "no cyclic recycling (declared)"),
    (r"^e0b\|setops\|port_error:call\|ValueError:parallel operations require equal lengths", DD,
     "e0b:README.md:103 (length 1 or n only); R recycles a length-1 y in GRanges pintersect", "no recycling (declared)"),
    (r"^e0b\|\*\|port_error:call\|rejects_names_padding", PB, "R pads short names with NA (not recycling)",
     "short names rejected instead of NA-padded"),
    (r"^e0b\|\*\|port_error:call\|rejects_double_arg", PB,
     "R truncates fractional coordinates/widths (IRanges(start=4, width=c(1.5,2.5,3.5)) -> widths 1 2 3); no e0b doc "
     "declares rejecting them", "fractional arguments rejected"),
    (r"^e0b\|(intra|construct|nearest|setops|inter|overlaps)\|(mismatch\|diff\|(int32_edge|maxgap_int32_max)|r_error_port_value\|int32_overflow)",
     DD, "e0b:README.md:52 'Integers: int64, no overflow check'. " + R_INT32, "int64 arithmetic (declared)"),
    (r"^e0b\|nearest\|mismatch\|diff\|(r_translation_dependent|coord_le0)$", RA,
     "probe_translate.py: R's answer changes when every range is shifted +1000 (2 cases on circular seqinfo not "
     "probed: zero-width subject at 0 / query at -2). " + R_SENTINEL, "R's GRanges sentinel collisions"),
    (r"^e0b\|nearest\|port_error:call\|ValueError:select must be one of \('first', 'all'\), got 'last'\|r_translation", PB,
     "follow(select='last') is R's API (precede: first/all, follow: last/all); e0b's follow accepts only first/all",
     "follow rejects select='last'"),
    (r"^e0b\|nearest\|r_error_port_value\|na_in_logic\|coord_le0$", RA, R_SENTINEL + "; R errors on valid input",
     "R fails on valid input with coordinates <= 0"),
    (r"^e0b\|nearest\|port_error:call\|ValueError:select must be one of \('first', 'all'\), got 'last'\|-$", PB,
     "R: follow(select=) is 'last' or 'all' (follow 'first' errors); e0b's follow takes 'first'/'all'",
     "follow uses precede's select vocabulary"),
    (r"^e0b\|nearest\|r_error_port_value\|other:'arg' should be one of", PB,
     "R: follow(select='first') errors (choices last/all); e0b accepts 'first'", "follow accepts select='first'"),
    (r"^e0b\|(nearest|overlaps|ordering|setops)\|port_error:call\|ValueError:(select|type|unknown ties_method 'a'|unknown resolve_empty 'start')", PB,
     "R match.arg partial matching ('a' -> 'all'/'any'/'average', 'start' -> 'start.x')", "no partial matching"),
    (r"^e0b\|ordering\|port_error:call\|ValueError:unknown ties_method 'last'", PB, "R rank(ties.method='last') is valid",
     "ties.method='last' unsupported"),
    (r"^e0b\|ordering\|r_error_port_value\|rank_default_ties", RA, R_RANK,
     "R's rank() default errors (S4 method defect); e0b has a working default (undeclared)"),
    (r"^e0b\|\*\|r_error_port_value\|(flag_validation|scalar_validation|type_validation)", PB,
     "R validates flags/scalars/types (isTRUEorFALSE, single integer); e0b coerces silently", "missing argument validation"),
    (r"^e0b\|(nearest|overlaps)\|r_error_port_value\|na_in_logic\|na_arg$", PB, "R errors on NA flags", "NA flag accepted"),
    (r"^e0b\|construct\|r_error_port_value\|na_input_rejected", PB, "R errors on NA seqnames / Seqinfo names", "accepts NA names"),
    (r"^e0b\|accessors\|r_error_port_value\|na_input_rejected", PB, "R: mcols(x) with NA names errors (rownames)", "no error"),
    (r"^e0b\|(overlaps|nearest|seqinfo)\|r_error_port_value\|(seqinfo_incompatible|other:The seqlevels)", PB,
     "R errors on incompatible seqinfo / invalid seqinfo<-", "accepts incompatible seqinfo"),
    (r"^e0b\|nearest\|port_error:call\|ValueError:conflicting", PB,
     "R precede/follow/distance do not merge seqinfo (no error; nearest does via findOverlaps); e0b raises",
     "raises on seqinfo conflicts R ignores in precede/follow/distance"),
    (r"^e0b\|setops\|port_error:call\|ValueError:parallel operations require (identical seqnames|compatible strands)", PB,
     "R GRanges pintersect flags non-hits (hit=FALSE), psetdiff passes x through for different seqnames/incompatible "
     "strands (e2:NOTES.md fuzz table); e0b NOTES marks its strand rule UNVERIFIED", "p* ops raise where R returns"),
    (r"^e0b\|setops\|port_error:call\|ValueError:some ranges do not overlap", PB,
     "R GRanges pintersect has no resolve.empty: a disjoint pair is hit=FALSE, zero width at x's start",
     "GRanges pintersect raises on disjoint pairs"),
    (r"^e0b\|setops\|r_error_port_value\|other:some ranges in 'y'", PB, "R psetdiff errors when y is strictly inside x",
     "no error"),
    (r"^e0b\|setops\|mismatch\|(mcols|mcols\+strand|strand|names)\|-$", PB,
     "R: GRanges pintersect adds a `hit` column, psetdiff/punion drop mcols, psetdiff takes y's strand for '*', "
     "ignore.strand keeps x's strand, IRanges p* use y's names when x has none; e0b NOTES call its p* metadata "
     "rules UNVERIFIED guesses", "p* metadata/strand rules"),
    (r"^e0b\|setops\|mismatch\|diff\|out_of_bounds$", PB, "R GRanges intersect/setdiff clip to [1, seqlength]", "no clipping"),
    (r"^e0b\|inter\|mismatch\|mcol_order\|-$", PB,
     "R revmap lists input indices in order(x) (start, width) within each reduced range; e0b NOTES: 'indices ascending "
     "... Order inside a revmap entry in R unknown. UNVERIFIED'", "revmap element order"),
    (r"^e0b\|inter\|mismatch\|mcols\|-$", PB,
     "R reduce(drop.empty.ranges=TRUE, with.revmap=TRUE) leaves dropped empty ranges out of revmap; e0b lists them",
     "revmap includes dropped empty ranges"),
    (r"^e0b\|inter\|mismatch\|diff\|(drop_empty_ranges|min_gapwidth_0)$", PB,
     "R: empty ranges bridge clusters under drop.empty.ranges=TRUE and identical empties merge at min.gapwidth=0 "
     "(e1:README.md 'R quirks reproduced'; e3:adr/0007, 0013)", "reduce empty-range rules"),
    (r"^e0b\|inter\|mismatch\|diff\|(zero_width|coord_le0|strand|multi_seqlevel|na_arg)$", PB,
     "R gaps: default window from non-empty ranges only, GRanges without seqlengths reports only strands with data, "
     "start/end vectors per seqlevel, end=NA means data extent; reduce(min.gapwidth=1.5) truncates to 1. e0b NOTES "
     "mark its gaps rule UNVERIFIED", "gaps/reduce bounds rules"),
    (r"^e0b\|inter\|port_error:call\|(IndexError|ValueError:operands|TypeError:int)", PB,
     "R gaps(start=/end=) take per-seqlevel vectors or NA; e0b crashes", "crash on vector/NA gap bounds"),
    (r"^e0b\|coverage\|mismatch\|diff\|na_arg$", PB,
     "R with NaN weights keeps NaN local (the port's actual input); e0b propagates NaN to later positions",
     "NaN weight poisons later positions"),
    (r"^e0b\|coverage\|mismatch\|diff\|zero_width$", PB,
     "R coverage default width counts zero-width ranges' ends (coverage(IRanges(c(1,20), width=c(3,0))) has length 19)",
     "zero-width ranges ignored for coverage length"),
    (r"^e0b\|coverage\|port_error:call\|ValueError:(weight must be numeric|The truth value)", PB,
     "R coverage(weight=<mcols column name>, width=<per-seqlevel vector>) work; e0b raises",
     "column-name weights / per-seqlevel widths unsupported (raises)"),
    (r"^e0b\|coverage\|port_error:call\|(IndexError|TypeError|UFuncTypeError|ValueError:operands)", PB,
     "R recycles/truncates shift and width vectors and accepts a shift column name; e0b crashes inside numpy",
     "coverage crashes on vector/double/column shift"),
    (r"^e0b\|coverage\|r_error_port_value\|other:", PB, "R rejects width=NA/TRUE and reordered named widths",
     "missing coverage validation"),
    (r"^e0b\|subset_concat\|port_error:call\|(mcols_type_conflict_crash|SchemaError)", PB,
     "R c() fills missing mcols columns with NA and coerces differing column types (int+dbl -> dbl, +chr -> chr); "
     "e0b refuses differing columns and polars raises SchemaError on differing types", "c() with differing mcols"),
    (r"^e0b\|subset_concat\|port_error:call\|AttributeError", PB, "R: c(<IRanges>, <IRanges>, 1L) coerces 1L", "crash"),
    (r"^e0b\|subset_concat\|r_error_port_value\|other:subscript contains invalid names", PB,
     "R: '' never matches a name in x[''] (error); e0b matches the empty name", "empty-string name subscript matches"),
    (r"^e0b\|construct\|port_error:call\|ValueError:at least two", PB,
     "R: IRanges(start=<int|lgl>) alone is the integer-set run form; IRanges() is empty", "start-only constructor raises"),
    (r"^e0b\|construct\|port_error:call\|ValueError:strand must be", PB, "R: NA strand becomes '*'", "NA strand rejected"),
    (r"^e0b\|construct\|port_error:call\|ValueError:seqlengths must be", PB, "R: Seqinfo(seqlengths=c(5, 6)) (double) is accepted",
     "double seqlengths rejected"),
    (r"^e0b\|construct\|port_error:call\|AttributeError", PB, "R: GRanges(seqinfo=<character>) is accepted", "crash"),
    (r"^e0b\|construct\|mismatch\|(mcols|names|seqlevel_order)\|-$", PB,
     "R keeps mcols of `ranges`, names whole doubles as '2' and logicals as 'TRUE', uses factor levels as seqlevels",
     "constructor metadata rules"),
    (r"^e0b\|accessors\|mismatch\|names\|-$", PB, "R mcols(x) has rownames = names(x)", "mcols() without rownames"),
    (r"^e0b\|accessors\|mismatch\|diff\|", PB,
     "R: mcols(<IRanges without columns>) is NULL, mcols(<GRanges without columns>) a 0-column DataFrame with "
     "nrow = length(x); e0b returns a 0x0 frame", "mcols() of an object without columns"),
    (r"^e0b\|overlaps\|mismatch\|names\|-$", PB, "R names countOverlaps by names(query)", "drops query names"),
    (r"^e0b\|overlaps\|mismatch\|diff\|(zero_width|plain|strand|coord_le0_r_invariant)$", PB,
     "probe_translate.py: R translation-invariant. R predicates: within with maxgap<=0 is unbounded, start/end/equal "
     "need no overlap (only |difference| <= maxgap), fractional minoverlap truncates, '-' strand end/start mirror; "
     "e0b NOTES mark its within/empty-range predicate UNVERIFIED", "overlap predicate differs from R"),
    (r"^e0b\|nearest\|mismatch\|(diff\|(zero_width|strand|plain|multi_seqlevel|coord_le0_r_invariant)|row_order\|-)$", PB,
     "probe_translate.py: R translation-invariant here. R: follow ties -> highest index (e3:adr/0010), '*' queries "
     "take the subject's strand direction, zero-width self queries, GRanges nearest(all) keeps only the precede side "
     "of a tie; e0b NOTES mark its strand rules UNVERIFIED", "nearest-family tie/strand rules differ from R"),
    (r"^e0b\|intra\|r_error_port_value\|other:solving row", PB, "R narrow(width=) alone errors", "accepts width-only narrow"),

    # ------------------------------------------------------------------ e2
    (r"^e2\|(overlaps|nearest|coverage)\|mismatch\|(equals_R_linear|diff)\|circular$", PB,
     R_CIRCULAR_PROBE + "; e2 documents no circular limitation (README Known limits lists only coordinates <= 0 "
     "and select ties)", "circular sequences treated as linear (undeclared)"),
    (r"^e2\|\*\|r_error_port_value\|length_validation\|circular$", PB,
     "R errors when coverage width exceeds a circular sequence's length; e2 has no circular handling (undeclared)",
     "circular coverage-width check missing"),
    (r"^e2\|coverage\|port_error:call\|(TypeError:weight must be numeric|ValueError:weight cannot contain NAs)\|na_weight_(int|dbl)$", DD,
     "e2:NOTES.md:27 'coverage: ... NA weight gives garbage (not pinned)'; the IRanges path refuses NA weights",
     "NA weights refused (R's NA-weight coverage judged garbage and deliberately not matched)"),
    (r"^e2\|coverage\|mismatch\|(nan_as_na\|na_weight_dbl|diff\|na_arg|diff\|strand)$", DD,
     "e2:NOTES.md:27 'NA weight gives garbage (not pinned)'; the GRanges path turns NA weights into NaN",
     "NA weights not matched to R (declared not pinned)"),
    (r"^e2\|nearest\|mismatch\|diff\|(r_translation_dependent|coord_le0)$", DD,
     "e2:README.md:109-112 'GRanges precede/follow/nearest with a range whose start or end is <= 0 ... can differ "
     "from R ... ranges_e2 implements the documented semantics'. " + R_SENTINEL, "geometric answer for coordinates <= 0 (declared)"),
    (r"^e2\|nearest\|r_error_port_value\|na_in_logic\|coord_le0$", DD, "e2:README.md:109-112 (coordinates <= 0)",
     "geometric answer where R errors (declared)"),
    (r"^e2\|nearest\|mismatch\|diff\|int32_edge$", RA, "probe/inputs: coordinates near -2^31 on one seqlevel; R's "
     "sentinel offsets overflow and return hits across seqnames (heldout/nearest/GRanges/117: chrM queries -> chrY subjects). "
     + R_SENTINEL, "R cross-sequence hits from int32 sentinel overflow"),
    (r"^e2\|nearest\|mismatch\|diff\|int32_edge_distance$", PB,
     "R distance() returns NA (with a warning) when the gap exceeds int32 (e1:NOTES.md 'distance() overflow -> NA'); "
     "e2 returns the int64 value, no doc", "distance overflow not NA"),
    (r"^e2\|nearest\|mismatch\|diff\|coord_le0_r_invariant$", PB,
     "probe_translate.py: R translation-invariant here (IRanges self nearest with zero-width ranges); e2 differs",
     "self-nearest with zero-width ranges"),
    (r"^e2\|ordering\|r_error_port_value\|rank_default_ties", DD, "e2:conformance/xfail.json 'ir.rank.default-is-error'; NOTES.md:104. "
     + R_RANK, "rank() default 'average' (declared)"),
    (r"^e2\|(setops|inter|overlaps)\|mismatch\|diff\|int32_edge$", RA, R_INT32, "R's int32 arithmetic gives a wrong answer"),
    (r"^e2\|\*\|r_error_port_value\|length_validation\|int32_edge$", RA, R_INT32, "R errors on valid input at 2^31-1"),
    (r"^e2\|intra\|(port_error:call\|int32_overflow|mismatch\|diff\|int32_edge)", RA, R_INT32_NA,
     "R returns NA coordinates its own validity check rejects; e2 raises / computes"),
    (r"^e2\|(intra|construct)\|r_error_port_value\|int32_overflow", PB,
     "R errors when a coordinate leaves int32 (IRanges constructor / IRanges shift); e2 builds the out-of-range object "
     "(its check_int32 is not applied there; no doc)", "no int32 bound check"),
    (r"^e2\|setops\|port_error:call\|resolve_empty_validation", RA, R_RESOLVE_EMPTY,
     "R validates resolve.empty only when an empty intersection occurs"),
    (r"^e2\|\*\|r_error_port_value\|(flag_validation|scalar_validation|type_validation)", PB,
     "R validates flags/scalars/types; e2 coerces (TRUE as width 1, NA flags as FALSE); no doc", "missing argument validation"),
    (r"^e2\|(intra|construct)\|r_error_port_value\|other:'(start|width|upstream|downstream|start' and 'width)", PB,
     "R requires numeric (not logical) start/width/upstream/downstream", "logical accepted as numbers"),
    (r"^e2\|intra\|r_error_port_value\|other:'(length\(shift\)|width' is longer|start', 'end' or 'width' is longer)", PB,
     "R errors when shift/width/narrow args are longer than x; e2 truncates", "over-long arguments accepted"),
    (r"^e2\|\*\|r_error_port_value\|length_validation\|-$", PB,
     "R errors when isCircular is shorter than seqnames (Seqinfo / isCircular<-); e2 recycles", "lenient Seqinfo lengths"),
    (r"^e2\|(nearest|overlaps)\|r_error_port_value\|na_in_logic\|na_arg$", PB, "R errors on NA ignore.strand/invert", "NA flag accepted"),
    (r"^e2\|\*\|port_error:call\|rejects_recycling", PB,
     "R recycles to zero length (IRanges(start=integer(0), width=1:3) is empty); e2 'cannot recycle a zero-length "
     "argument', no doc", "zero-length recycling rejected"),
    (r"^e2\|construct\|port_error:call\|ValueError:seqnames must have the length", PB,
     "R recycles a length-1 ranges over seqnames (GRanges(c('chr1','chr2'), IRanges(1, width=5)))", "length-1 ranges not recycled"),
    (r"^e2\|construct\|port_error:call\|ValueError:at least two", PB, "R: IRanges(start=<int|lgl>) alone is the run form",
     "start-only constructor raises"),
    (r"^e2\|construct\|port_error:call\|AttributeError", PB, "R: GRanges(seqinfo=<character>) is accepted", "crash"),
    (r"^e2\|construct\|mismatch\|names\|-$", PB, "R names whole doubles '2', logicals 'TRUE'", "str() of names"),
    (r"^e2\|intra\|port_error:call\|ValueError:(start|end|width) cannot contain NAs", PB,
     "R narrow(start=NA)/end=NA/width=NA means 'not given'; e2 rejects NA", "NA narrow arguments rejected"),
    (r"^e2\|nearest\|port_error:call\|ValueError:sequence chr\w+ has incompatible", PB,
     "R precede/follow/distance do not merge seqinfo and return values; e2 raises", "seqinfo merge check R does not do here"),
    (r"^e2\|seqinfo\|port_error:call\|ValueError:seqlevels to be dropped", PB,
     "R seqlevels<- with a same-length value renames positionally (chr2 -> chr2_r); e2 NOTES 'rename if nothing in "
     "common' misses partial renames", "seqlevels<- partial rename"),
    (r"^e2\|seqinfo\|r_error_port_value\|other:The seqlevels", PB, "R refuses a reordered Seqinfo in seqinfo<-", "accepted"),
    (r"^e2\|(nearest|overlaps|ordering)\|port_error:call\|ValueError:(select|type|ties_method) must be", PB,
     "R match.arg partial matching ('a' -> 'all'/'any'/'average')", "no partial matching"),
    (r"^e2\|overlaps\|mismatch\|diff\|plain$", PB, "R truncates fractional minoverlap (5.5 -> 5)", "fractional minoverlap not truncated"),
    (r"^e2\|overlaps\|mismatch\|names\|-$", PB, "R names countOverlaps by names(query); no e2 doc", "drops query names"),
    (r"^e2\|subset_concat\|mismatch\|mcols\|-$", PB, "R as.character(double) uses 15 significant digits", "double->chr formatting"),
    (r"^e2\|subset_concat\|mismatch\|mcol_types", PB, "R: column only in an empty object is logical NA", "takes the empty object's type"),
    (r"^e2\|subset_concat\|port_error:call\|TypeError", PB, "R: c(<IRanges>, <IRanges>, 1L) coerces 1L", "c() rejects integer"),
    (r"^e2\|subset_concat\|r_error_port_value\|other:subscript contains invalid names", PB, "R: '' never matches a name",
     "empty-string name subscript matches"),
    (r"^e2\|setops\|mismatch\|(seqinfo|names|mcols)\|-$", PB,
     "R: GRanges pintersect keeps x's seqlevels only and replaces an existing `hit` column; IRanges p* take y's names "
     "when x has none", "p* metadata rules"),
    (r"^e2\|accessors\|mismatch\|names\|-$", PB, "R mcols(x) has rownames = names(x); no e2 doc", "mcols() without rownames"),
    (r"^e2\|accessors\|mismatch\|diff\|", PB, "R: mcols(<IRanges without columns>) is NULL; e2 returns an empty frame",
     "mcols() of an IRanges without columns"),
    (r"^e2\|accessors\|r_error_port_value\|na_input_rejected", PB, "R mcols(x) errors when names(x) has NA", "no error"),
    (r"^e2\|inter\|mismatch\|mcols\|-$", PB, "R reduce(min.gapwidth=0) drops repeated empties from revmap (e1:NOTES 00:20)",
     "revmap keeps repeated empties"),
    (r"^e2\|inter\|mismatch\|diff\|zero_width$", PB, "R truncates min.gapwidth=1.5 to 1", "fractional min.gapwidth"),
    (r"^e2\|coverage\|port_error:call\|(TypeError:'float'|TypeError:shift must|TypeError:weight must be numeric\|-)", PB,
     "R accepts coverage(width=1.5) (truncated) and shift/weight as an mcols column name", "coverage argument forms"),
    (r"^e2\|coverage\|r_error_port_value\|other:", PB, "R rejects width=NA/TRUE, over-long shift, reordered named width",
     "missing coverage validation"),

    # ------------------------------------------------------------------ e3
    (r"^e3\|(overlaps|nearest|coverage)\|mismatch\|(equals_R_linear|diff)\|circular$", DD,
     "e3:adr/0014-circular-metadata-only.md; README.md 'Circular sequences get no wrap-around (ADR-0014)'. " + R_CIRCULAR_PROBE,
     "no circular wrap-around (declared)"),
    (r"^e3\|\*\|r_error_port_value\|length_validation\|circular$", DD, "e3:adr/0014-circular-metadata-only.md",
     "circular coverage-width check not implemented (declared)"),
    (r"^e3\|\*\|port_error:call\|rejects_recycling", DD,
     "e3:adr/0006-refuse-bad-lengths.md; spec/10-iranges-core.md:13 CORE.3 'other length mismatches are refused (R warns "
     "and recycles; ADR-0006)'", "no cyclic recycling (declared)"),
    (r"^e3\|.*\|port_error:call\|rejects_scalar_genome", PB,
     "R recycles a scalar genome over seqnames (Seqinfo(c('a','b','c'), genome='hg')); e3's own CORE.3/ADR-0006 say "
     "length-1 values recycle", "length-1 genome rejected (contradicts its own rule)"),
    (r"^e3\|\*\|port_error:call\|rejects_names_padding", PB, "R pads short names with NA (not recycling)", "short names rejected"),
    (r"^e3\|(intra|construct|setops|inter|overlaps)\|(mismatch\|diff\|(int32_edge|maxgap_int32_max)|r_error_port_value\|int32_overflow)", DD,
     "e3:adr/0005-int64.md 'Results beyond int32 range are valid in Python and are not refused'. " + R_INT32,
     "int64 arithmetic (declared)"),
    (r"^e3\|\*\|r_error_port_value\|length_validation\|int32_edge$", DD, "e3:adr/0005-int64.md. " + R_INT32, "int64 (declared)"),
    (r"^e3\|nearest\|mismatch\|diff\|int32_edge_distance$", DD, "e3:adr/0005-int64.md (R returns NA when the distance "
     "exceeds int32)", "int64 distance (declared)"),
    (r"^e3\|nearest\|mismatch\|diff\|(r_translation_dependent|int32_edge|coord_le0)$", RA,
     "probe_translate.py: R's answer moves with the coordinates (2 cases on circular seqinfo unprobed; zero-width subject "
     "at 0 / query at -2). e3 declares only the zero-width-at-1 corner (spec/95-granges.md:60). " + R_SENTINEL,
     "R's GRanges sentinel collisions"),
    (r"^e3\|nearest\|port_error:call\|ValueError:select must be 'first' or 'all'\|r_translation", PB,
     "R follow(select=) is 'last'/'all'; e3 spec/70-nearest.md NEAR.4 uses 'first' for follow", "follow select vocabulary"),
    (r"^e3\|nearest\|r_error_port_value\|na_in_logic\|coord_le0$", RA, R_SENTINEL + "; R errors on valid input",
     "R fails on valid input with coordinates <= 0"),
    (r"^e3\|nearest\|mismatch\|diff\|(multi_seqlevel|strand|zero_width)$", DD,
     "e3:adr/0011-granges-nearest-ties.md (GRanges nearest(all) keeps all minimal-distance ties; R keeps only the precede "
     "side); spec/70-nearest.md NEAR.2 (zero-width self queries implementation-defined)", "GRanges nearest ties (declared)"),
    (r"^e3\|nearest\|mismatch\|diff\|coord_le0_r_invariant$", PB,
     "probe_translate.py: R translation-invariant. Self nearest(x): NEAR.2 declares only zero-width self queries "
     "implementation-defined, but positive-width queries differ too (rnd0-nearest-self q2=[2,3]: R 36, e3 5)",
     "self-nearest candidate set differs beyond the declared zero-width corner"),
    (r"^e3\|nearest\|port_error:call\|ValueError:select must be 'first' or 'all'\|-$", PB,
     "R: follow(select=) is 'last'|'all' (follow 'first' errors \"'arg' should be one of last, all\"); e3 "
     "spec/70-nearest.md NEAR.4 and adr/0010 assume R's follow takes 'first'", "follow uses the wrong select vocabulary"),
    (r"^e3\|nearest\|r_error_port_value\|other:'arg' should be one of", PB, "R errors on follow(select='first')",
     "accepts follow(select='first')"),
    (r"^e3\|(nearest|overlaps|ordering)\|port_error:call\|ValueError:(select must be one of|type must be one of|ties_method must be 'average')", PB,
     "R match.arg partial matching ('a') and rank(ties.method='last') is valid R", "choice arguments: no partial matching / no 'last'"),
    (r"^e3\|ordering\|r_error_port_value\|rank_default_ties", RA, R_RANK, "R's rank() default errors; e3 defaults to average"),
    (r"^e3\|setops\|port_error:call\|(resolve_empty_validation|ValueError:The truth value)", RA, R_RESOLVE_EMPTY,
     "R validates resolve.empty only when needed"),
    (r"^e3\|\*\|r_error_port_value\|(flag_validation|scalar_validation|type_validation)", PB,
     "R validates flags/scalars/types; e3 coerces", "missing argument validation"),
    (r"^e3\|(intra|construct|coverage)\|r_error_port_value\|other:", PB,
     "R requires logical-without-NA flank(start=), numeric start/width/upstream/downstream, valid coverage width",
     "missing argument validation"),
    (r"^e3\|(nearest|overlaps)\|r_error_port_value\|na_in_logic\|na_arg$", PB, "R errors on NA flags", "NA flag accepted"),
    (r"^e3\|construct\|r_error_port_value\|na_input_rejected", PB,
     "R errors on NA/empty seqnames; e3 stores NA as the string 'None'", "NA seqnames accepted as 'None'"),
    (r"^e3\|accessors\|r_error_port_value\|na_input_rejected", PB, "R mcols(x) errors when names(x) has NA", "no error"),
    (r"^e3\|\*\|port_error:call\|na_arg_crash", PB,
     "R gaps(start=NA)/end=NA means 'from the data'; e3 receives a masked array (adapters/e3.md 'A length-1 NA argument "
     "is passed as a masked array') and crashes with MaskError", "NA gap bounds crash"),
    (r"^e3\|intra\|port_error:call\|ValueError:(start|end|width) contains NA", PB, "R narrow(start=NA) etc. means 'not given'",
     "NA narrow arguments rejected"),
    (r"^e3\|\*\|port_error:call\|rejects_double_arg", PB,
     "R truncates fractional coordinates/shifts (IRanges(start=c(1.5,2.5,3.5), width=1) -> 1 2 3); e3 CORE.1 says "
     "'integer sequences' but declares no divergence", "fractional arguments rejected"),
    (r"^e3\|construct\|port_error:call\|ValueError:exactly two", PB,
     "R accepts a consistent start/end/width triple and the start-only run form; e3 spec CORE.1 claims 'exactly two' "
     "as R behaviour (no ADR)", "constructor argument sets narrower than R"),
    (r"^e3\|construct\|port_error:call\|ValueError:strand values", PB, "R: NA strand becomes '*'", "NA strand rejected"),
    (r"^e3\|construct\|mismatch\|(mcols|names|seqlevel_order)\|-$", PB,
     "R keeps mcols(ranges) (e3 spec GR.1 says 'its mcols are ignored', no ADR), names whole doubles '2' / 'TRUE', "
     "uses factor levels as seqlevels", "constructor metadata rules"),
    (r"^e3\|accessors\|mismatch\|names\|-$", PB, "R mcols(x) rownames = names(x); ranges(<GRanges>) keeps names(x)",
     "accessors drop names"),
    (r"^e3\|accessors\|mismatch\|diff\|", PB, "R: mcols(<IRanges without columns>) is NULL, GRanges a 0-column frame with nrow",
     "mcols() of an object without columns (adapters/e3.md notes the 0x0 frame)"),
    (r"^e3\|overlaps\|mismatch\|names\|-$", PB, "R names countOverlaps by names(query)", "drops query names"),
    (r"^e3\|setops\|mismatch\|(mcols|seqinfo\+mcols|names)\|-$", PB,
     "R: GRanges pintersect adds a `hit` column and keeps x's seqlevels; IRanges p* take y's names when x has none",
     "p* metadata rules"),
    (r"^e3\|setops\|mismatch\|diff\|out_of_bounds$", PB,
     "R clips GRanges intersect/setdiff to [1, seqlength] and drops positions < 1 even with unknown length; e3 spec "
     "GR.12 clips only when the length is known", "no clipping below 1 / unknown lengths"),
    (r"^e3\|subset_concat\|mismatch\|(mcols|mcol_types)", PB,
     "R formats doubles with 15 significant digits in chr columns; a column only in an empty object is logical NA",
     "c() mcols coercion"),
    (r"^e3\|subset_concat\|r_error_port_value\|other:subscript contains invalid names", PB, "R: '' never matches a name",
     "empty-string name subscript matches"),
    (r"^e3\|inter\|port_error:call\|TypeError:Only length", PB, "R gaps(start=/end=) take per-seqlevel vectors",
     "crash on vector gap bounds"),
    (r"^e3\|coverage\|port_error:call\|ValueError:(The truth value|setting an array)", PB,
     "R coverage(width=<per-seqlevel vector | NA>) works; e3 crashes", "coverage width forms crash"),
    (r"^e3\|nearest\|port_error:call\|ValueError:sequence 'chr\w+' has incompatible", PB,
     "R precede/follow/distance do not merge seqinfo; e3 raises", "seqinfo merge check R does not do here"),

    # ------------------------------------------------------------------ e4
    (r"^e4\|(overlaps|nearest)\|mismatch\|(equals_R_linear|diff)\|circular$", PB,
     R_CIRCULAR_PROBE + "; e4 implements circular wrap in coverage only (README.md:67, NOTES.md:16) and documents no "
     "linear treatment of circular sequences in overlaps/nearest", "circular overlaps/nearest computed linearly (undeclared)"),
    (r"^e4\|(coverage|inter)\|mismatch\|diff\|circular$", PB,
     "probe_circular.py: circularity matters in R here; e4's circular coverage/gaps disagree with R", "circular corner"),
    (r"^e4\|coverage\|mismatch\|diff\|(na_arg|coord_le0|strand)$", PB,
     "e4 receives NA weights as masked arrays (adapters/e4.md 'RVector with NA -> masked array') and counts them as 0; "
     "R's (local-NA) answer here is NA (coverage(IRanges(c(1,2,5,8,11),width=c(1,3,3,3,2)), weight=c(NA,3,1,NA,4)) "
     "= NA 3 1 NA 4); no e4 doc", "NA weights treated as 0"),
    (r"^e4\|\*\|port_error:call\|rejects_recycling", DD,
     "e4:adr/0007-vector-args-recycle-from-one.md 'length-1 recycles; any other length != n raises' (zero length too)",
     "no cyclic recycling (declared)"),
    (r"^e4\|\*\|port_error:call\|rejects_names_padding", PB, "R pads short names with NA (not recycling)", "short names rejected"),
    (r"^e4\|(intra|construct|setops|inter|overlaps|nearest)\|(mismatch\|diff\|(int32_edge|maxgap_int32_max)|r_error_port_value\|int32_overflow|port_error:call\|int32_overflow)", DD,
     "e4:adr/0006-int32-overflow-raises.md, adr/0010-r-overflow-defects-not-ported.md. " + R_INT32,
     "int32 overflow raises / R overflow defects not ported (declared)"),
    (r"^e4\|\*\|r_error_port_value\|length_validation\|int32_edge$", DD, "e4:adr/0010-r-overflow-defects-not-ported.md "
     "(ir-edge-disjoin-end-int-max)", "R overflow defect not ported (declared)"),
    (r"^e4\|nearest\|mismatch\|diff\|int32_edge_distance$", DD, "e4:adr/0010 (ir-edge-distance-overflow-is-na pinned)", "declared"),
    (r"^e4\|nearest\|mismatch\|diff\|(r_translation_dependent|coord_le0)$", RA,
     "probe_translate.py: R's answer moves with the coordinates (circular-seqinfo cases unprobed). " + R_SENTINEL,
     "R's GRanges sentinel collisions"),
    (r"^e4\|nearest\|r_error_port_value\|na_in_logic\|coord_le0$", RA, R_SENTINEL + "; R errors on valid input",
     "R fails on valid input with coordinates <= 0"),
    (r"^e4\|construct\|port_error:call\|ValueError:at least two", DD,
     "e4:adr/0011-iranges-needs-two-of-three.md 'fewer than two of start/end/width is always a ValueError'",
     "start-only / empty IRanges() refused (declared)"),
    (r"^e4\|overlaps\|mismatch\|names\|-$", DD, "e4:adr/0009-count-overlaps-drops-names.md", "countOverlaps unnamed (declared)"),
    (r"^e4\|ordering\|r_error_port_value\|rank_default_ties", RA, R_RANK, "R's rank() default errors"),
    (r"^e4\|setops\|port_error:call\|(resolve_empty_validation|ValueError:The truth value)", RA, R_RESOLVE_EMPTY,
     "R validates resolve.empty only when needed"),
    (r"^e4\|\*\|port_error:call\|rejects_double_arg", PB,
     "R truncates fractional coordinates/shifts/widths; e4 raises 'non-integer value' (no ADR)", "fractional arguments rejected"),
    (r"^e4\|\*\|r_error_port_value\|(flag_validation|scalar_validation|type_validation)", PB, "R validates flags/scalars/types",
     "missing argument validation"),
    (r"^e4\|(intra|construct|coverage)\|r_error_port_value\|other:", PB,
     "R requires logical-without-NA flank(start=), numeric start/width/upstream/downstream, valid coverage width/shift",
     "missing argument validation"),
    (r"^e4\|\*\|r_error_port_value\|length_validation\|-$", PB, "R errors on short Seqinfo seqlengths/isCircular", "recycled"),
    (r"^e4\|(nearest|overlaps)\|r_error_port_value\|na_in_logic\|na_arg$", PB, "R errors on NA flags", "NA flag accepted"),
    (r"^e4\|accessors\|r_error_port_value\|na_input_rejected", PB, "R mcols(x) errors when names(x) has NA", "no error"),
    (r"^e4\|(intra|coverage)\|port_error:call\|ValueError:(start|end|width): NA not allowed", PB,
     "R narrow(start=NA)/coverage(width=NA) mean 'not given'", "NA arguments rejected"),
    (r"^e4\|construct\|port_error:call\|ValueError:strand values", PB, "R: NA strand becomes '*'", "NA strand rejected"),
    (r"^e4\|construct\|port_error:call\|ValueError:argument lengths differ", DD, "e4:adr/0007-vector-args-recycle-from-one.md",
     "constructor arguments recycle only from length 1 (declared)"),
    (r"^e4\|construct\|mismatch\|(mcols|names|seqlevel_order)\|-$", PB,
     "R keeps mcols(ranges), names whole doubles '2' / 'TRUE', uses factor levels as seqlevels", "constructor metadata"),
    (r"^e4\|accessors\|mismatch\|names\|-$", PB, "R mcols(x) rownames = names(x)", "mcols() without rownames"),
    (r"^e4\|accessors\|mismatch\|diff\|", PB, "R: mcols(<IRanges without columns>) is NULL; GRanges 0-col frame with nrow",
     "mcols() of an object without columns (adapters/e4.md notes the 0x0 frame)"),
    (r"^e4\|setops\|mismatch\|(seqinfo|strand|names)\|-$", PB,
     "R: pintersect keeps x's seqlevels; psetdiff '*' takes y's strand (e2:NOTES GR.PSETDIFF.STRAND); IRanges p* take "
     "y's names when x has none", "p* metadata/strand rules"),
    (r"^e4\|setops\|mismatch\|diff\|out_of_bounds$", PB, "R clips GRanges intersect/setdiff to [1, seqlength]", "no clipping"),
    (r"^e4\|inter\|mismatch\|diff\|(zero_width|coord_le0|na_arg)$", PB,
     "R disjoin splits ranges at zero-width ranges (disjoin-revmap: 1-10 split at 5), gaps(start=1.5) truncates, "
     "gaps(end=NA) means the data extent; e4:NOTES.md 'disjoin ignores zero-width' is a misread of R", "disjoin/gaps rules"),
    (r"^e4\|inter\|mismatch\|diff\|(drop_empty_ranges|min_gapwidth_0)$", PB,
     "R: empties bridge under drop.empty.ranges=TRUE; identical empties merge at min.gapwidth=0 (e3:adr/0007, 0013)",
     "reduce empty-range rules"),
    (r"^e4\|inter\|mismatch\|mcols\|-$", PB, "R reduce/disjoin(with.revmap=TRUE) on an empty GRanges has an empty revmap column",
     "missing revmap column on empty input"),
    (r"^e4\|inter\|port_error:call\|TypeError", PB, "R gaps(start=/end=) take per-seqlevel vectors", "vector gap bounds rejected"),
    (r"^e4\|subset_concat\|mismatch\|(mcols|mcol_types)", PB, "R formats doubles with 15 significant digits; empty-only column is lgl",
     "c() mcols coercion"),
    (r"^e4\|subset_concat\|port_error:call\|AttributeError", PB, "R: c(<IRanges>, <IRanges>, 1L) coerces 1L", "crash"),
    (r"^e4\|nearest\|port_error:call\|ValueError:sequence chr\w+ has incompatible", PB,
     "R precede/follow/distance do not merge seqinfo", "seqinfo merge check R does not do here"),
    (r"^e4\|(nearest|overlaps|ordering)\|port_error:call\|ValueError:'(select|type|ties_method)' must be one of", PB,
     "R match.arg partial matching ('a')", "no partial matching"),
    (r"^e4\|nearest\|mismatch\|diff\|(multi_seqlevel|strand)$", PB,
     "R GRanges nearest(select='all') keeps only the precede side of an exact tie (e2:NOTES GR.NEAREST.TIE)", "nearest tie rule"),
    (r"^e4\|coverage\|port_error:call\|(TypeError|ValueError)", PB,
     "R accepts coverage(shift=1.5) (truncated), weight/shift as mcols column names", "coverage argument forms"),

    # ------------------------------------------------------------------ biocpy (iranges 0.7.3, genomicranges 0.9.0)
    (r"^biocpy\|overlaps\|mismatch\|select_result_as_hits", AF,
     "adapters/biocpy.py:312 `return self._selected_hits(res, n_query) or _Hits(...)`: for an empty query the "
     "selection is RVector('int', []), falsy through RVector.__len__ (runner/adapter_base.py:63), so the projectable "
     "integer(0) falls back to Hits (e4/iranges-overlaps/ir-fo-select-empty: R integer(0), port Hits 0x5)",
     "empty select= result projected as Hits"),
    (r"^biocpy\|overlaps\|mismatch\|select_not_applied", PB,
     "BiocPy find_overlaps(select='first'|'last'|'arbitrary') returns several subjects for one query (e1/grid_circular/"
     "circ-long0-fo-first: query 7 -> subjects 1 and 3; R first 1, last 3); the adapter cannot project it "
     "(adapters/biocpy.py `_selected_hits` returns None) and keeps Hits", "select= not applied by the port"),
    (r"^biocpy\|overlaps\|mismatch\|hits_not_sorted_by_query", PB,
     "rule R2 (spec/COMPARATOR.md) requires Hits in non-decreasing query order; BiocPy find_overlaps returns them "
     "sorted by subject (heldout/findOverlaps/IRanges/58: R from 2 3 4 4 5 5, port 2 3 4 5 4 5); same pairs",
     "findOverlaps Hits not sorted by query"),
    (r"^biocpy\|overlaps\|mismatch\|diff\|int32_edge$", PB,
     "heldout/findOverlaps/GRanges/63 (maxgap=0, select='last'): R and e0-e4 give subject 2 for queries 6-9 "
     "(chrY:2147483615-2147483623 vs chrY:2147483616-2147483618 '*'); BiocPy only for query 8",
     "overlaps near 2^31 missed"),
    (r"^biocpy\|nearest\|mismatch\|diff\|int32_edge_distance$", PB,
     "R ?distance: NA when the distance exceeds int32; BiocPy returns 0 for ranges ~2^31 apart "
     "(e4/iranges-edge/ir-edge-distance-overflow-is-na)", "distance beyond int32 is 0"),
    (r"^biocpy\|intra\|mismatch\|diff\|port_int32_wrap$", PB,
     "e1/granges_intra/gr-shift-overflow-error: shift(x, 2147483647) gives BiocPy starts -2147483639... (int32 "
     "wrap-around); R's NA starts are invalid too (" + R_INT32_NA + "), but a silently wrapped coordinate is the port's own "
     "defect", "int32 wrap-around to negative coordinates"),
    (r"^biocpy\|(overlaps|nearest|coverage|setops|seqinfo|subset_concat|inter|ordering)\|mismatch\|equals_R_linear\|circular$", PB,
     R_CIRCULAR_PROBE + "; BiocPy documents no circular limitation", "circular sequences treated as linear"),
    (r"^biocpy\|nearest\|(mismatch\|diff\|(r_translation_dependent|int32_edge)|port_error:call\|IndexError:only integers.*\|r_translation_dependent)", RA,
     "probe_translate.py: R's answer moves with the coordinates / overflows near 2^31. " + R_SENTINEL,
     "R's GRanges sentinel collisions"),
    (r"^biocpy\|nearest\|r_error_port_value\|na_in_logic\|coord_le0$", RA, R_SENTINEL, "R fails on valid input"),
    (r"^biocpy\|(setops|inter)\|mismatch\|diff\|int32_edge$", RA, R_INT32, "R's int32 arithmetic gives a wrong answer"),
    (r"^biocpy\|intra\|mismatch\|diff\|int32_edge$", RA, R_INT32_NA + "; BiocPy gives the int64 arithmetic answer "
     "(promoters(IRanges(2147483642, width=3), 0, 100): width 100)", "R's overflow result is not a contract"),
    (r"^biocpy\|\*\|r_error_port_value\|length_validation\|int32_edge$", RA, R_INT32, "R errors on valid input at 2^31-1"),
    (r"^biocpy\|setops\|port_error:call\|(resolve_empty_validation|ValueError:The truth value)", RA, R_RESOLVE_EMPTY, "lazy R validation"),
    (r"^biocpy\|coverage\|mismatch\|num_type", PB,
     "BiocPy coverage returns float64 arrays for integer coverage (R integer Rle); adapter projects dtype faithfully",
     "coverage is double"),
    (r"^biocpy\|overlaps\|mismatch\|(num_type|num_type\+names)\|-$", PB,
     "BiocPy count_overlaps returns float64 (R integer); names differ", "countOverlaps type"),
    (r"^biocpy\|nearest\|port_error:call\|IndexError:only integers", PB,
     "GenomicRanges.precede/follow crash inside the port (genomicranges/GenomicRanges.py:2578, 2658 `ret_result[q]`); "
     "reproduced directly in the port venv", "GRanges precede/follow crash"),
    (r"^biocpy\|\*\|port_error:call\|rejects_recycling", PB,
     "R recycles start/width/names/mcols (and NA-pads names); BiocPy requires equal lengths, undocumented", "no recycling"),
    (r"^biocpy\|(seqinfo|setops|subset_concat|inter|overlaps|ordering)\|mismatch\|(seqinfo_fields|seqinfo|seqinfo\+mcols)", PB,
     "BiocPy merge_SeqInfo / set ops / c() drop seqlengths, isCircular and genome that R merges (seqinfo-merge-basic: "
     "R lengths 1,5,7; port 1,NA,7) and return an empty Seqinfo for empty results", "Seqinfo merge loses fields"),
    (r"^biocpy\|(seqinfo|setops|subset_concat|inter|overlaps|nearest)\|r_error_port_value\|seqinfo_incompatible", PB,
     "R errors on incompatible seqlengths/circularity/genome; BiocPy merges silently (keeps NA)", "no seqinfo conflict check"),
    (r"^biocpy\|coverage\|port_error:call\|TypeError:int\(\) argument", PB,
     "every member is IRanges coverage (array weights are its documented API; scalar-only GRanges weight vectors are "
     "Unsupported since adapter fef1051); the adapter passes NA as None and the port crashes converting None to int. "
     "9 members are R-order-sensitive (probe_coverage_na order_dependent) but R's value at the given order is the "
     "local-NA sum, so the port is held to it", "NA weights crash"),
    (r"^biocpy\|(coverage|intra)\|port_error:call\|Exception:input length N exceeds expected length N", PB,
     "heldout/coverage/IRanges/41: coverage(<empty IRanges>, shift=4) is integer(0) in R; BiocPy's array-valued "
     "IRanges.coverage normalizes the scalar shift to length 1 and rejects it against 0 ranges (also width=69 in "
     "heldout/coverage/IRanges/17, resize of an empty GRanges in heldout/resize/GRanges/43)",
     "scalar argument on empty input raises"),
    (r"^biocpy\|accessors\|r_error_port_value\|na_input_rejected", PB,
     "e1/construct/ir-mcols-accessor-na-names: R mcols(IRanges(names=c('a',NA))) errors (NA rownames); BiocPy has "
     "already turned the NA name into the string 'None' and returns the frame. mcols(x, use.names=FALSE) succeeds in R "
     "but asks for different semantics, so it is no equivalent-call evidence of an R accident",
     "NA name stored as 'None'; no NA-rownames error"),
    (r"^biocpy\|coverage\|mismatch\|diff\|na_arg$", PB, "BiocPy counts masked (NA) weights as 0; R's local answer is NA",
     "NA weights treated as 0"),
    (r"^biocpy\|", PB,
     "BiocPy's READMEs promise the Bioconductor semantics ('Python implementation of the IRanges Bioconductor package'; "
     "inclusive ends 'to be consistent with Bioconductor') and declare no divergence; the examined cases show the "
     "port's own behaviour (validation, crashes, metadata, strand/zero-width rules) differing from R",
     "BiocPy behaviour differs from R (no declared reason)"),
]
