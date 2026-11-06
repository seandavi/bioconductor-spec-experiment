# Eval oracle: run cases in R (normative) and emit canonical projections.
#   Rscript oracle/oracle.R <in.jsonl> <out.jsonl>
# in : one {id, inputs, call, universe_call} per line (spec/CASES.md; universe_call may be null)
# out: one record per line, flushed per case:
#   {id, status: ok|excluded, reason?, inputs, result, r_class, sorted_by_query, universe?, warnings, domain}
# Contract: spec/CASES.md (decode), spec/PROJECTION.md (shapes). Uses jsonlite only to PARSE.
suppressPackageStartupMessages(suppressWarnings({
  library(IRanges); library(GenomicRanges); library(jsonlite)
}))
options(warn = 1)

# ---- JSON writer (exact doubles: jsonlite digits=NA is only 15 digits) -------------------
jstr <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  s <- gsub("[\001-\037]", "", s)
  paste0("\"", s, "\"")
}
jnum <- function(x) {
  if (is.na(x)) return(if (is.nan(x)) "\"NaN\"" else "null")
  if (is.infinite(x)) return(if (x > 0) "\"Inf\"" else "\"-Inf\"")
  s <- sprintf("%.17g", x)
  if (grepl("^-?[0-9]+$", s)) s <- paste0(s, ".0")
  s
}
js <- function(x) {
  if (is.null(x)) return("null")
  if (is.list(x)) {
    if (length(x) == 0) return(if (is.null(names(x))) "[]" else "{}")
    items <- vapply(x, js, "")
    if (is.null(names(x))) return(paste0("[", paste(items, collapse = ","), "]"))
    return(paste0("{", paste(vapply(names(x), jstr, ""), items, sep = ":", collapse = ","), "}"))
  }
  if (length(x) != 1) stop("js: non-scalar atomic; wrap with arr()")
  if (is.character(x)) return(if (is.na(x)) "null" else jstr(x))
  if (is.logical(x)) return(if (is.na(x)) "null" else if (x) "true" else "false")
  if (is.integer(x)) return(if (is.na(x)) "null" else as.character(x))
  if (is.double(x)) return(jnum(x))
  stop("js: unsupported scalar ", class(x)[1])
}
arr <- function(x) { x <- unname(x); lapply(seq_along(x), function(i) x[[i]]) }
obj <- function(...) list(...)

# ---- eval decode (spec/CASES.md "Literals") -----------------------------------------------
NAv <- function(e) if (is.null(e)) NA else e
vec_decode <- function(v, type) {
  vals <- lapply(v[[type]], NAv)
  x <- switch(type,
    int = as.integer(unlist(vals)), dbl = as.numeric(unlist(vals)),
    chr = as.character(unlist(vals)), lgl = as.logical(unlist(vals)))
  if (length(vals) == 0) x <- switch(type, int = integer(), dbl = numeric(), chr = character(), lgl = logical())
  if (!is.null(v$names)) names(x) <- as.character(unlist(lapply(v$names, NAv)))
  x
}
is_vec_lit <- function(v) is.list(v) && !is.null(names(v)) && any(names(v) %in% c("int", "dbl", "chr", "lgl"))
is_obj_lit <- function(v, k) is.list(v) && length(v) == 1 && identical(names(v), k)

decode <- function(v, env) {
  if (is.null(v)) return(NULL)
  if (is_vec_lit(v)) return(vec_decode(v, intersect(names(v), c("int", "dbl", "chr", "lgl"))[1]))
  if (is.list(v) && length(v) == 1 && identical(names(v), "$ref")) {
    nm <- v[["$ref"]]
    if (!exists(nm, envir = env$inputs, inherits = FALSE)) stop("unknown input: ", nm)
    return(get(nm, envir = env$inputs))
  }
  if (is.list(v) && length(v) == 1 && identical(names(v), "$call")) return(run_call(v[["$call"]], env)$value)
  if (is.list(v) && !is.null(names(v))) {
    k <- names(v)[1]; b <- v[[1]]
    g <- function(f) if (is.null(b[[f]])) NULL else decode(b[[f]], env)
    return(switch(k,
      IRanges = { x <- IRanges(start = g("start"), end = g("end"), width = g("width"), names = g("names"))
                  with_mcols(x, b$mcols, env) },
      GRanges = { ir <- IRanges(start = g("start"), end = g("end"), width = g("width"), names = g("names"))
                  a <- list(seqnames = g("seqnames"), ranges = ir)
                  if (!is.null(b$strand)) a$strand <- g("strand")
                  if (!is.null(b$seqinfo)) a$seqinfo <- if (!is.null(b$seqinfo[["$ref"]])) decode(b$seqinfo, env) else seqinfo_decode(b$seqinfo, env)
                  with_mcols(do.call(GRanges, a), b$mcols, env) },
      Seqinfo = seqinfo_decode(b, env),
      DataFrame = {
        cols <- columns_decode(b$columns, env)
        if (!is.null(b$nrow)) { df <- S4Vectors::make_zero_col_DFrame(b$nrow); for (nm in names(cols)) df[[nm]] <- cols[[nm]]; df }
        else do.call(DataFrame, c(cols, list(check.names = FALSE)))
      },
      Factor = { lv <- as.character(unlist(b$levels)); cd <- as.integer(unlist(lapply(b$codes, NAv)))
                 factor(lv[cd], levels = lv) },
      IntegerList = { l <- lapply(b$elements, function(e) as.integer(unlist(lapply(e, NAv))))
                      if (!is.null(b$names)) names(l) <- as.character(unlist(lapply(b$names, NAv)))
                      IntegerList(l) },
      stop("unknown literal: ", k)))
  }
  # bare scalars: chr / lgl / int / dbl of length 1
  v
}
seqinfo_decode <- function(b, env) {
  g <- function(f) if (is.null(b[[f]])) NULL else decode(b[[f]], env)
  a <- list(seqnames = g("seqnames"))
  for (f in c("seqlengths", "isCircular", "genome")) if (!is.null(b[[f]])) a[[f]] <- g(f)
  do.call(Seqinfo, a)
}
columns_decode <- function(cols, env) {
  out <- lapply(cols, function(p) decode(p[[2]], env))
  names(out) <- vapply(cols, function(p) p[[1]], "")
  out
}
with_mcols <- function(x, mc, env) {
  if (is.null(mc) || length(mc) == 0) return(x)
  cols <- columns_decode(mc, env)
  mcols(x) <- do.call(DataFrame, c(cols, list(check.names = FALSE)))
  x
}

# ---- evaluation ---------------------------------------------------------------------------
run_call <- function(cl, env) {
  args <- lapply(cl$args, decode, env = env)
  kwargs <- lapply(cl$kwargs, decode, env = env)
  first <- if (length(args)) args[1] else if (length(kwargs)) kwargs[1] else list(NULL)
  env$domain <- domain_of(first[[1]])  # recorded even if the call below errors
  list(value = do.call(cl$op, c(args, kwargs)), first = first[[1]])
}

# ---- projection (spec/PROJECTION.md) ---------------------------------------------------------
names_of <- function(x) { n <- names(x); if (is.null(n) || length(x) == 0) NULL else arr(n) }
Unproj <- function(cls) stop(structure(class = c("unprojectable", "error", "condition"),
  list(message = paste0("r_unprojectable:", cls), call = NULL)))

p_vec <- function(x, with_names = TRUE) {
  type <- if (is.integer(x)) "int" else if (is.double(x)) "dbl" else if (is.character(x)) "chr" else if (is.logical(x)) "lgl" else Unproj(class(x)[1])
  obj(kind = "vector", type = type, values = arr(x), names = if (with_names) names_of(x) else NULL)
}
p_factor <- function(x) obj(kind = "Factor", levels = arr(levels(x)), codes = arr(as.integer(x)))
p_intlist <- function(x) obj(kind = "IntegerList", elements = unname(lapply(as.list(x), function(e) arr(as.integer(e)))), names = names_of(x))
p_col <- function(col) {
  if (is.factor(col)) return(p_factor(col))
  if (is(col, "IntegerList")) return(p_intlist(col))
  if (is.atomic(col) && !is.object(col)) return(p_vec(col, with_names = FALSE))
  Unproj(class(col)[1])
}
p_mcols <- function(x) {
  m <- mcols(x, use.names = FALSE)
  if (is.null(m) || ncol(m) == 0) return(list())
  lapply(seq_len(ncol(m)), function(j) list(colnames(m)[j], p_col(m[[j]])))
}
p_seqinfo <- function(si) obj(kind = "Seqinfo", seqnames = arr(seqlevels(si)), seqlengths = arr(seqlengths(si)),
  isCircular = arr(isCircular(si)), genome = arr(genome(si)))
p_rle_body <- function(x) {
  rv <- runValue(x)
  type <- if (is.integer(rv)) "int" else if (is.double(rv)) "dbl" else if (is.logical(rv)) "lgl" else Unproj(paste0("Rle<", class(rv)[1], ">"))
  list(type = type, values = arr(rv), lengths = arr(runLength(x)))
}

project <- function(x) {
  if (is.null(x)) return(obj(kind = "NULL"))
  if (is(x, "GenomicRanges")) {
    return(obj(kind = "GRanges", seqnames = arr(as.character(seqnames(x))), start = arr(start(x)), width = arr(width(x)),
      strand = arr(as.character(strand(x))), names = names_of(x), mcols = p_mcols(x), seqinfo = p_seqinfo(seqinfo(x))))
  }
  if (is(x, "IntegerRanges")) {
    return(obj(kind = "IRanges", start = arr(start(x)), width = arr(width(x)), names = names_of(x), mcols = p_mcols(x)))
  }
  if (is(x, "Hits")) {
    return(obj(kind = "Hits", from = arr(queryHits(x)), to = arr(subjectHits(x)), nLnode = nLnode(x), nRnode = nRnode(x)))
  }
  if (is(x, "Seqinfo")) return(p_seqinfo(x))
  if (is(x, "RleList")) {
    return(obj(kind = "RleList", names = names_of(x), elements = unname(lapply(as.list(x), p_rle_body))))
  }
  if (is(x, "Rle")) {
    if (is.factor(runValue(x))) return(p_vec(as.character(x)))
    b <- p_rle_body(x)
    return(obj(kind = "Rle", type = b$type, values = b$values, lengths = b$lengths))
  }
  if (is(x, "IntegerList")) return(p_intlist(x))
  if (is(x, "DataFrame")) {
    cols <- lapply(seq_len(ncol(x)), function(j) list(colnames(x)[j], p_col(x[[j]])))
    return(obj(kind = "DataFrame", nrow = nrow(x), rownames = names_of(structure(seq_len(nrow(x)), names = rownames(x))), columns = cols))
  }
  if (is.factor(x)) return(p_factor(x))
  if (is.atomic(x) && !is.object(x)) return(p_vec(x))
  Unproj(class(x)[1])
}

domain_of <- function(v) {
  if (is(v, "GenomicRanges")) "GRanges" else if (is(v, "IntegerRanges")) "IRanges" else if (is(v, "Seqinfo")) "Seqinfo" else "other"
}

# ---- one case -------------------------------------------------------------------------------
run_case <- function(case) {
  warns <- character()
  env <- new.env(); env$inputs <- new.env(); env$domain <- "other"
  rec <- list(id = case$id)
  setTimeLimit(elapsed = 120, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  out <- withCallingHandlers(
    tryCatch({
      # inputs: an R failure here excludes the case
      in_proj <- tryCatch({
        for (nm in names(case$inputs)) assign(nm, decode(case$inputs[[nm]], env), envir = env$inputs)
        p <- lapply(names(case$inputs), function(nm) project(get(nm, envir = env$inputs)))
        names(p) <- names(case$inputs)
        p
      }, error = function(e) structure(list(msg = conditionMessage(e)), class = "input_error"))
      if (inherits(in_proj, "input_error")) {
        list(status = "excluded", reason = "r_input_error", detail = substr(in_proj$msg, 1, 300))
      } else {
        res <- tryCatch({
          rc <- run_call(case$call, env)
          list(ok = TRUE, value = rc$value, domain = domain_of(rc$first))
        }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
        warn_snapshot <- warns
        if (!res$ok) {
          list(status = "ok", inputs = in_proj, result = obj(kind = "error", message = res$message),
               r_class = NULL, domain = env$domain)
        } else {
          proj <- tryCatch(project(res$value), unprojectable = function(e) e)
          if (inherits(proj, "unprojectable")) {
            list(status = "excluded", reason = conditionMessage(proj), detail = class(res$value)[1])
          } else {
            r <- list(status = "ok", inputs = in_proj, result = proj, r_class = class(res$value)[1], domain = env$domain)
            if (is(res$value, "Hits")) r$sorted_by_query <- is(res$value, "SortedByQueryHits")
            if (!is.null(case$universe_call) && proj$kind == "vector") {
              u <- tryCatch(project(run_call(case$universe_call, env)$value), error = function(e) NULL)
              if (!is.null(u)) r$universe <- u
            }
            r
          }
        }
      }
    }, error = function(e) list(status = "excluded", reason = "r_harness_error", detail = substr(conditionMessage(e), 1, 300))),
    warning = function(w) warns <<- c(warns, conditionMessage(w)))
  out$warnings <- arr(unique(warns))
  c(rec, out)
}

args <- commandArgs(trailingOnly = TRUE)
lines <- readLines(args[1], warn = FALSE)
con <- file(args[2], "w")
for (ln in lines) {
  if (!nzchar(ln)) next
  case <- fromJSON(ln, simplifyVector = FALSE)
  rec <- tryCatch(run_case(case), error = function(e)
    list(id = case$id, status = "excluded", reason = "r_harness_error", detail = substr(conditionMessage(e), 1, 300), warnings = list()))
  writeLines(js(rec), con)
  flush(con)
}
close(con)
