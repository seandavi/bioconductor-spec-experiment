# Usage: Rscript run_cases.R in.json out.json
# in.json : [{id, inputs, call | raw}, ...]   (see tests/differential/README.md)
# out.json: {versions, results: {id: {value|error, warnings}}}
suppressWarnings(suppressMessages({
  library(IRanges); library(GenomicRanges); library(jsonlite)
}))
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
source(file.path(here, "encode.R"))
files <- commandArgs(trailingOnly = TRUE)

# ---- literal decoding -------------------------------------------------------
vec_lit <- function(v, type = NULL) {
  if (length(v) == 0) {
    return(switch(if (is.null(type)) "int" else type,
      int = integer(), num = numeric(), chr = character(), lgl = logical()))
  }
  x <- unlist(lapply(v, function(e) if (is.null(e)) NA else e))
  switch(if (is.null(type)) "" else type,
    int = as.integer(x), num = as.numeric(x), chr = as.character(x),
    lgl = as.logical(x), x)
}

df_lit <- function(cols) {
  cols <- lapply(cols, decode)
  do.call(DataFrame, c(cols, list(check.names = FALSE)))
}

seqinfo_lit <- function(l) {
  g <- function(k) if (is.null(l[[k]])) NULL else decode(l[[k]])
  a <- list(seqnames = g("seqnames"))
  if (!is.null(l$seqlengths)) a$seqlengths <- as.integer(g("seqlengths"))
  if (!is.null(l$is_circular)) a$isCircular <- g("is_circular")
  if (!is.null(l$genome)) a$genome <- g("genome")
  do.call(Seqinfo, a)
}

iranges_lit <- function(l) {
  g <- function(k) if (is.null(l[[k]])) NULL else decode(l[[k]])
  x <- IRanges(start = g("start"), end = g("end"), width = g("width"), names = g("names"))
  if (!is.null(l$mcols)) mcols(x) <- df_lit(l$mcols)
  x
}

granges_lit <- function(l) {
  g <- function(k) if (is.null(l[[k]])) NULL else decode(l[[k]])
  ir <- IRanges(start = g("start"), end = g("end"), width = g("width"), names = g("names"))
  si <- if (is.null(l$seqinfo)) NULL else seqinfo_lit(l$seqinfo)
  st <- g("strand"); if (is.null(st)) st <- "*"
  x <- GRanges(seqnames = g("seqnames"), ranges = ir, strand = st, seqinfo = si)
  if (!is.null(l$mcols)) mcols(x) <- df_lit(l$mcols)
  x
}

decode <- function(v) {
  if (is.null(v)) return(NA)
  if (is.list(v) && !is.null(names(v))) {
    k <- names(v)[1]
    return(switch(k,
      IRanges = iranges_lit(v[[1]]), GRanges = granges_lit(v[[1]]),
      Seqinfo = seqinfo_lit(v[[1]]), DataFrame = df_lit(v[[1]]),
      int = vec_lit(v[[1]], "int"), num = vec_lit(v[[1]], "num"),
      chr = vec_lit(v[[1]], "chr"), lgl = vec_lit(v[[1]], "lgl"),
      stop("unknown literal: ", k)))
  }
  if (is.list(v)) return(vec_lit(v))
  v
}

# ---- evaluation -------------------------------------------------------------
resolve <- function(v, env) {
  if (is.character(v) && length(v) == 1 && startsWith(v, "$")) {
    if (startsWith(v, "$$")) return(substring(v, 2))
    nm <- substring(v, 2)
    if (!exists(nm, envir = env, inherits = FALSE)) stop("unknown input: ", nm)
    return(get(nm, envir = env))
  }
  if (is.list(v) && !is.null(names(v)) && identical(names(v), "call")) {
    return(run_call(v$call, env))
  }
  decode(v)
}

run_call <- function(cl, env) {
  a <- lapply(cl$args, resolve, env = env)
  k <- lapply(cl$kwargs, resolve, env = env)
  do.call(cl$fn, c(a, k))
}

run_case <- function(case) {
  warns <- character()
  env <- new.env()
  out <- withCallingHandlers(
    tryCatch({
      for (nm in names(case$inputs)) assign(nm, decode(case$inputs[[nm]]), envir = env)
      val <- if (!is.null(case$raw)) eval(parse(text = case$raw$r), envir = env)
             else run_call(case$call, env)
      list(value = enc(val))
    }, error = function(e) list(error = u(conditionMessage(e)))),
    warning = function(w) {
      # record only: muffling changes control flow in IRanges (warning -> error paths)
      warns <<- c(warns, conditionMessage(w))
    })
  out$warnings <- warns
  out
}

cases <- fromJSON(readLines(files[1], warn = FALSE) |> paste(collapse = "\n"),
  simplifyVector = FALSE)
res <- list()
for (case in cases) res[[case$id]] <- run_case(case)
ver <- paste0("R ", getRversion(), "; IRanges ", packageVersion("IRanges"),
  "; GenomicRanges ", packageVersion("GenomicRanges"),
  "; S4Vectors ", packageVersion("S4Vectors"))
writeLines(toJSON(list(versions = u(ver), results = res), null = "null", na = "null",
  auto_unbox = FALSE, digits = NA, pretty = FALSE), files[2])
