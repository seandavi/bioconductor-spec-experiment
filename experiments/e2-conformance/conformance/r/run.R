# R runner.
#   Rscript conformance/r/run.R validate [group...]
#   Rscript conformance/r/run.R generate [group...]   # writes expected/<group>.json (normative)
#   Rscript conformance/r/run.R check [group...]      # re-runs R vs expected; writes reports/r.json
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
root <- dirname(here)  # conformance/
source(file.path(here, "codec.R"))

# ---- loading + validation -------------------------------------------------

load_groups <- function(groups = character()) {
  files <- sort(list.files(file.path(root, "cases"), pattern = "\\.json$", full.names = TRUE))
  all <- lapply(files, read_json)
  names(all) <- vapply(all, function(g) g$group, "")
  if (length(groups)) {
    missing <- setdiff(groups, names(all))
    if (length(missing)) stop("unknown group(s): ", paste(missing, collapse = ", "))
    all <- all[groups]
  }
  all
}

find_refs <- function(node) {
  if (!is.list(node)) return(character())
  own <- if (!is.null(names(node)) && !is.null(node[["$ref"]])) node[["$ref"]] else character()
  c(own, unlist(lapply(node, find_refs)))
}

validate_groups <- function(all_groups) {
  errs <- character()
  ids <- character()
  for (g in all_groups) {
    clause_used <- setNames(rep(FALSE, length(g$clauses)), names(g$clauses))
    for (cs in g$cases) {
      if (cs$id %in% ids) errs <- c(errs, paste("duplicate id", cs$id))
      ids <- c(ids, cs$id)
      if (is.null(cs$clause) || !(cs$clause %in% names(g$clauses)))
        errs <- c(errs, paste(cs$id, ": unknown clause", cs$clause))
      else clause_used[[cs$clause]] <- TRUE
      if (is.null(cs$op) || !nzchar(cs$op)) errs <- c(errs, paste(cs$id, ": missing op"))
      for (r in find_refs(list(cs$args, cs$kwargs)))
        if (!(r %in% names(g$fixtures))) errs <- c(errs, paste(cs$id, ": unresolved ref", r))
    }
    for (cl in names(clause_used)[!clause_used])
      errs <- c(errs, paste(g$group, ": clause without case", cl))
  }
  errs
}

# ---- running --------------------------------------------------------------

decode_fixtures <- function(g) {
  fx <- list()
  for (nm in names(g$fixtures)) fx[[nm]] <- decode(g$fixtures[[nm]], fx)
  fx
}

run_case <- function(cs, fx) {
  args <- lapply(cs$args, decode, fx = fx)
  kwargs <- lapply(cs$kwargs, decode, fx = fx)
  warns <- character()
  res <- tryCatch(
    withCallingHandlers(
      list(value = do.call(cs$op, c(args, kwargs))),
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
      message = function(m) invokeRestart("muffleMessage")),
    error = function(e) list(error = list(message = conditionMessage(e), class = as.list(class(e)))))
  if (!is.null(res$error)) return(list(error = res$error))
  list(value = encode(res$value), warnings = as.list(warns))
}

normalize <- function(tree) from_json(to_json(tree))

meta <- function() {
  pv <- function(p) as.character(packageVersion(p))
  list(R = R.version.string, IRanges = pv("IRanges"), GenomicRanges = pv("GenomicRanges"),
       S4Vectors = pv("S4Vectors"), Seqinfo = pv("Seqinfo"),
       generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
}

# ---- commands -------------------------------------------------------------

cmd_generate <- function(all_groups) {
  for (g in all_groups) {
    fx <- decode_fixtures(g)
    out <- list()
    for (cs in g$cases) out[[cs$id]] <- run_case(cs, fx)
    doc <- list(meta = meta(), cases = out)
    writeLines(to_json(doc), file.path(root, "expected", paste0(g$group, ".json")))
    n_err <- sum(vapply(out, function(o) !is.null(o$error), TRUE))
    n_warn <- sum(vapply(out, function(o) length(o$warnings) > 0, TRUE))
    cat(sprintf("%-18s %3d cases (%d errors, %d with warnings)\n",
                g$group, length(out), n_err, n_warn))
  }
}

cmd_check <- function(all_groups) {
  report <- list()
  bad <- 0
  for (g in all_groups) {
    path <- file.path(root, "expected", paste0(g$group, ".json"))
    if (!file.exists(path)) stop("no expected file for ", g$group, "; run generate")
    expected <- read_json(path)$cases
    fx <- decode_fixtures(g)
    for (cs in g$cases) {
      exp <- expected[[cs$id]]
      st <- tryCatch({
        got <- normalize(run_case(cs, fx))
        if (is.null(exp)) "fail"
        else if (!is.null(exp$error)) { if (!is.null(got$error)) "pass" else "fail" }
        else if (is.null(got$error) && tree_equal(got$value, exp$value)) {
          # codec self-consistency on the stored golden
          rt <- if (is.null(exp$value)) NULL else normalize(encode(decode(exp$value)))
          if (tree_equal(rt, exp$value)) "pass" else "fail"
        } else "fail"
      }, error = function(e) "error")
      if (st != "pass") { bad <- bad + 1; cat("R", st, cs$id, "\n") }
      report[[cs$id]] <- st
    }
  }
  if (length(commandArgs(TRUE)) <= 1)  # full run only
    writeLines(to_json(list(implementation = "R", cases = report)),
               file.path(root, "reports", "r.json"))
  cat(sprintf("R check: %d cases, %d not passing\n", length(report), bad))
  quit(status = as.integer(bad > 0))
}

main <- function() {
  a <- commandArgs(TRUE)
  cmd <- if (length(a)) a[1] else "check"
  groups <- load_groups(a[-1])
  errs <- validate_groups(groups)
  if (length(errs)) { cat(paste("INVALID:", errs), sep = "\n"); quit(status = 2) }
  switch(cmd, validate = cat("valid:", sum(vapply(groups, function(g) length(g$cases), 1L)), "cases\n"),
         generate = cmd_generate(groups), check = cmd_check(groups), stop("unknown command ", cmd))
}
main()
