#!/usr/bin/env Rscript
# Usage: Rscript runners/r/run_cases.R [--check]
# Reads cases/<domain>/<group>.yaml, writes goldens/<domain>/<group>.json.
# --check recomputes and exits 1 if committed goldens differ (drift detector).
# R is normative: nothing in a golden is hand-edited.

suppressPackageStartupMessages({
  library(IRanges); library(GenomicRanges); library(yaml); library(jsonlite)
})
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
root <- normalizePath(file.path(here, "..", ".."))
source(file.path(here, "build.R"))
source(file.path(here, "project.R"))
check <- "--check" %in% commandArgs(TRUE)

run_call <- function(call, inputs, kwargs_override = NULL) {
  args <- unname(lapply(call$args, function(n) inputs[[n]]))
  kw <- lapply(call$kwargs, resolve_kwarg, inputs = inputs)
  for (n in names(kwargs_override)) kw[[n]] <- kwargs_override[[n]]
  do.call(call$op, c(args, kw))
}

# Per-element allowed values for implementation-defined `select` results,
# computed by R from the select="all" Hits.
candidates_from <- function(h, n) {
  lapply(seq_len(n), function(i) {
    s <- subjectHits(h)[queryHits(h) == i]
    if (length(s)) s else NA_integer_
  })
}

run_case <- function(case) {
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers({
      inputs <- build_inputs(case$inputs)
      out <- project(run_call(case$call, inputs))
      if (!is.null(case$oneof_from)) {
        h <- run_call(case$call, inputs, case$oneof_from$kwargs)
        out$candidates <- candidates_from(h, queryLength(h))
      }
      out
    }, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) list(kind = unbox("error"), message = unbox(conditionMessage(e))))
  res$warnings <- as.character(warns)
  res
}

versions <- list(R = as.character(getRversion()),
  IRanges = as.character(packageVersion("IRanges")),
  GenomicRanges = as.character(packageVersion("GenomicRanges")),
  S4Vectors = as.character(packageVersion("S4Vectors")))
header <- lapply(versions, unbox)

to_json <- function(x) toJSON(x, auto_unbox = FALSE, na = "null", digits = NA, pretty = TRUE)

files <- sort(Sys.glob(file.path(root, "cases", "*", "*.yaml")))
seen <- character(0)
drift <- FALSE
for (f in files) {
  domain <- basename(dirname(f)); group <- sub("\\.yaml$", "", basename(f))
  cases <- yaml.load_file(f)$cases
  out <- list()
  for (case in cases) {
    if (case$id %in% seen) stop("duplicate case id: ", case$id)
    seen <- c(seen, case$id)
    out[[case$id]] <- run_case(case)
  }
  gold <- file.path(root, "goldens", domain, paste0(group, ".json"))
  new <- to_json(list(generated_by = header, cases = out))
  if (check) {
    old <- if (file.exists(gold)) paste(readLines(gold, warn = FALSE), collapse = "\n") else ""
    if (!identical(as.character(new), old)) {
      cat("DRIFT:", gold, "\n"); drift <- TRUE
    }
  } else {
    dir.create(dirname(gold), showWarnings = FALSE, recursive = TRUE)
    writeLines(new, gold)
    cat("wrote", length(out), "cases ->", sub(paste0(root, "/"), "", gold), "\n")
  }
}
if (check) { if (drift) quit(status = 1) else cat("goldens up to date\n") }
