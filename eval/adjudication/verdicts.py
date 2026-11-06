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
R_EMPTY_TYPE = ("GRanges precede/follow return SortedByQueryHits when the query or subject is empty (integer "
                "otherwise; IRanges: integer either way) and GRanges nearest returns logical NA for an empty "
                "subject (IRanges: integer NA): R's return type depends on emptiness")
R_RESOLVE_EMPTY = ("pintersect(x, x, resolve.empty=NA|TRUE|'a'|1:2|integer(0)) returns; the same values error once an "
                   "empty intersection must be resolved (pintersect(x, y) with a non-overlapping pair)")
R_CIRCULAR_PROBE = ("probe_circular.py: with every isCircular set FALSE R returns exactly the port's value, i.e. the "
                    "port computes linear semantics on a circular sequence")
HF_NAN = ("runner/adapter_base.py:82,88 RVector.np() maps a dbl NA to NaN (only int/lgl/chr raise; ADAPTERS.md "
          "'dbl NA and NaN both become nan') and adapters use it for weights (adapters/e0b.py:76; e0.py:64-67 mirrors "
          "it), while spec/PROJECTION.md and COMPARATOR exact equality score NaN != NA. The port's answer equals R's "
          "answer for the NaN input it was given: coverage(IRanges(c(1,2,5,8,11),width=c(1,3,3,3,2)), "
          "weight=c(NaN,3,1,NaN,4)) is NaN 3 1 NaN 4")
R_TRIM_NAMES = "trim(x, use.names=FALSE) keeps names(x) (R ignores the argument); shift(use.names=FALSE) drops them"

RULES: list[tuple[str, str, str, str]] = [
    # ------------------------------------------------------------------ cross-port R accidents
    (r"\|coverage\|.*\|r_nonlocal_na$", RA, R_COV_NA,
     "R's NA-weight coverage is not a sum here (NA leaks past the NA range / depends on order)"),
    (r"\|nearest\|mismatch\|diff\|empty$", RA, R_EMPTY_TYPE,
     "empty query/subject: R's GRanges result class (Hits / logical) differs from its own non-empty and IRanges answers"),
    (r"^(e0|e0b|e2|e3|e4|biocpy)\|intra\|mismatch\|diff\|strand$", RA, R_RESIZE_FIX,
     "R's ignore.strand=TRUE resize ignores the per-range fix vector"),
    (r"^(e0|e0b|e2|e3|e4|biocpy)\|intra\|port_error:call\|ValueError:'fix' must", RA, R_RESIZE_EMPTY,
     "R skips fix validation on empty input only"),
    # ------------------------------------------------------------------ e1
    (r"^e1\|coverage\|.*\|r_nonlocal_na$", RA, R_COV_NA,
     "R's NA-weight coverage is not a sum here; e1 reproduces the defect but cannot be held to it"),
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
    (r"^e1\|(overlaps|setops|inter)\|mismatch\|diff\|int32_edge$", DD, "e1:adrs/0005-coordinate-overflow-is-an-error.md; " + R_INT32,
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
    (r"^(e0|e0b)\|coverage\|mismatch\|nan_as_na\|na_weight_dbl$", HF, HF_NAN,
     "double NA reaches the port as NaN by the harness's own conversion, then NaN is scored against NA"),
    (r"^e0\|coverage\|mismatch\|(nan_as_na\+)?num_type\|na_weight_int$", AF,
     "adapters/e0.py:64-67 converts an integer vector holding NA to a float64 NaN array; e0's as_weights "
     "(src/ranges_e0/_coverage.py:25-28) then computes double coverage. R got an integer weight",
     "integer NA weights reach the port as doubles: the type change is the adapter's"),
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
    (r"^e0\|coverage\|port_error:call\|mcols_column_name_arg", AF,
     "R's coverage(weight=\"col\") names an mcols column; e0's coverage takes numbers only (README API table, "
     "_coverage.py:13 as_weights); adapters/e0.py passes the string through instead of raising Unsupported",
     "adapter passes a column-name weight the port's API does not have"),
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
    (r"^e0\|\*\|port_error:call\|na_arg_crash", AF,
     "adapters/e0.py:64-67 passes R NA as NaN, which the e0 port documents no meaning for; gaps then fails "
     "converting NaN to int. R: gaps(x, start=NA) means 'use the data'", "NA bound mistranslated to NaN"),
    (r"^e0\|inter\|mismatch\|diff\|na_arg$", AF,
     "adapters/e0.py:14 maps NA in lgl to None, which is e0's 'argument absent' default (seqlength bound); "
     "R's gaps(end=NA) means 'max end of the data', a different call", "NA bound mistranslated to the default"),
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
    (r"^e0\|(setops|inter|overlaps|nearest|intra)\|mismatch\|diff\|int32_edge$", DD,
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
    (r"^e0b\|(intra|construct|nearest|setops|inter|overlaps)\|(mismatch\|diff\|int32_edge|r_error_port_value\|int32_overflow)",
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
]
