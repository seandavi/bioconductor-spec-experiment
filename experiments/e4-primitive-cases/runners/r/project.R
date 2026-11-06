# Project R results onto primitive fields (see PLAN.md 4.2). The Python
# runner mirrors this function. Scalars are unboxed; vectors stay arrays.

null_json <- function() unbox(NA)  # serialises as JSON null
nz <- function(x) if (is.null(x)) null_json() else x
empty_obj <- function() structure(list(), names = character(0))

proj_vec <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.null(x)) return(list())
  as.vector(x)
}

proj_mcols <- function(x) {
  mc <- mcols(x)
  out <- empty_obj()
  if (is.null(mc)) return(out)
  for (n in colnames(mc)) {
    col <- mc[[n]]
    out[[n]] <- if (is(col, "List")) lapply(as.list(col), as.vector) else proj_vec(col)
  }
  out
}

proj_seqinfo <- function(si) {
  list(seqnames = seqnames(si), seqlengths = proj_vec(seqlengths(si)),
       isCircular = proj_vec(isCircular(si)), genome = proj_vec(genome(si)))
}

project <- function(x) {
  if (is(x, "GRanges")) {
    list(kind = unbox("GRanges"),
         seqnames = as.character(seqnames(x)), start = start(x), end = end(x),
         width = width(x), strand = as.character(strand(x)),
         names = nz(names(x)), mcols = proj_mcols(x),
         seqinfo = proj_seqinfo(seqinfo(x)))
  } else if (is(x, "IntegerRanges")) {
    list(kind = unbox("IRanges"),
         start = start(x), end = end(x), width = width(x),
         names = nz(names(x)), mcols = proj_mcols(x))
  } else if (is(x, "Hits")) {
    list(kind = unbox("Hits"), query = queryHits(x), subject = subjectHits(x),
         query_length = unbox(queryLength(x)), subject_length = unbox(subjectLength(x)))
  } else if (is(x, "Rle")) {
    list(kind = unbox("Rle"), values = proj_vec(runValue(x)), lengths = runLength(x))
  } else if (is(x, "RleList")) {
    list(kind = unbox("RleList"), names = nz(names(x)),
         elements = unname(lapply(as.list(x), function(r)
           list(values = proj_vec(runValue(r)), lengths = runLength(r)))))
  } else if (is(x, "Seqinfo")) {
    c(list(kind = unbox("Seqinfo")), proj_seqinfo(x))
  } else if (is.atomic(x)) {
    list(kind = unbox("vector"), type = unbox(typeof(x)),
         values = proj_vec(x), names = nz(names(x)))
  } else stop("cannot project class ", class(x)[1])
}
