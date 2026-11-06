# Build native R objects from primitive input specs (see PLAN.md 4.1).
# YAML quirks handled here: `null` inside a sequence arrives as a list with
# NULL; `[]` arrives as list(); both normalised by `vec()`.

vec <- function(x, type = NULL) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    x <- unlist(lapply(x, function(e) if (is.null(e)) NA else e))
    if (is.null(x)) x <- logical(0)
  }
  switch(if (is.null(type)) "none" else type,
    int = as.integer(x), num = as.numeric(x), chr = as.character(x),
    lgl = as.logical(x), x)
}

build_mcols <- function(x, mc) {
  if (is.null(mc)) return(x)
  cols <- lapply(mc, vec)
  mcols(x) <- S4Vectors::DataFrame(cols, check.names = FALSE)
  x
}

build_iranges <- function(s) {
  args <- list()
  for (k in c("start", "end", "width")) if (!is.null(s[[k]])) args[[k]] <- vec(s[[k]], "int")
  if (!is.null(s$names)) args$names <- vec(s$names, "chr")
  build_mcols(do.call(IRanges, args), s$mcols)
}

build_seqinfo <- function(s) {
  args <- list(seqnames = vec(s$seqnames, "chr"))
  if (!is.null(s$seqlengths)) args$seqlengths <- vec(s$seqlengths, "int")
  if (!is.null(s$isCircular)) args$isCircular <- vec(s$isCircular, "lgl")
  if (!is.null(s$genome)) args$genome <- vec(s$genome, "chr")
  do.call(Seqinfo, args)
}

build_granges <- function(s, inputs) {
  ir <- build_iranges(s[intersect(names(s), c("start", "end", "width", "names"))])
  args <- list(seqnames = vec(s$seqnames, "chr"), ranges = ir)
  if (!is.null(s$strand)) args$strand <- vec(s$strand, "chr")
  si <- s$seqinfo
  if (is.character(si)) args$seqinfo <- inputs[[si]]
  else if (!is.null(si)) args$seqinfo <- build_seqinfo(si)
  build_mcols(do.call(GRanges, args), s$mcols)
}

build_input <- function(spec, inputs) {
  stopifnot(length(spec) == 1)
  kind <- names(spec)
  s <- spec[[1]]
  switch(kind,
    IRanges = build_iranges(s),
    GRanges = build_granges(s, inputs),
    Seqinfo = build_seqinfo(s),
    int = vec(s, "int"), num = vec(s, "num"),
    chr = vec(s, "chr"), lgl = vec(s, "lgl"),
    stop("unknown input kind: ", kind))
}

build_inputs <- function(specs) {
  inputs <- list()
  for (n in names(specs)) inputs[[n]] <- build_input(specs[[n]], inputs)
  inputs
}

# kwargs: scalars pass through, sequences become typed vectors,
# {ref: name} passes a built input.
resolve_kwarg <- function(v, inputs) {
  if (is.list(v) && identical(names(v), "ref")) return(inputs[[v$ref]])
  if (is.list(v)) return(vec(v))
  v
}
