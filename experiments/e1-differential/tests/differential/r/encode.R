# R value -> canonical JSON-ready list. Keep in sync with adapter.py (Python side).
# Conventions: indices stay 1-based (R's); vectors are always arrays;
# factors / factor-Rles are expanded to character vectors.

u <- jsonlite::unbox

vec_type <- function(x) {
  if (is.logical(x)) "lgl"
  else if (is.integer(x)) "int"
  else if (is.numeric(x)) "num"
  else if (is.character(x)) "chr"
  else stop("unsupported vector type: ", class(x)[1])
}

enc_vec <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  nm <- names(x)
  list(class = u("vector"), type = u(vec_type(x)),
       values = unname(as.vector(x)),
       names = if (is.null(nm)) NULL else unname(nm))
}

is_atomic_vec <- function(x) is.atomic(x) && !is.null(x)

enc_mcols <- function(df) {
  if (is.null(df) || ncol(df) == 0L) return(NULL)
  cols <- lapply(seq_len(ncol(df)), function(i) enc(df[[i]]))
  names(cols) <- colnames(df)
  cols
}

enc_names <- function(x) {
  nm <- names(x)
  if (is.null(nm)) NULL else unname(nm)
}

enc_seqinfo <- function(si) {
  list(class = u("Seqinfo"),
       seqnames = seqlevels(si),
       seqlengths = unname(seqlengths(si)),
       is_circular = unname(isCircular(si)),
       genome = unname(genome(si)))
}

enc_list_elem <- function(e) if (is_atomic_vec(e) && !is.factor(e)) unname(as.vector(e)) else enc(e)

enc <- function(x) {
  if (is.null(x)) return(NULL)
  if (is(x, "GRanges")) {
    return(list(class = u("GRanges"),
      seqnames = as.character(seqnames(x)),
      start = start(x), width = width(x),
      strand = as.character(strand(x)),
      names = enc_names(x),
      mcols = enc_mcols(mcols(x)),
      seqinfo = enc_seqinfo(seqinfo(x))))
  }
  if (is(x, "IRanges")) {
    return(list(class = u("IRanges"), start = start(x), width = width(x),
      names = enc_names(x), mcols = enc_mcols(mcols(x))))
  }
  if (is(x, "Seqinfo")) return(enc_seqinfo(x))
  if (is(x, "Hits")) {
    return(list(class = u("Hits"), from = from(x), to = to(x),
      nLnode = u(nLnode(x)), nRnode = u(nRnode(x))))
  }
  if (is(x, "RleList")) {
    return(list(class = u("RleList"), names = enc_names(x),
      elements = lapply(as.list(x), enc)))
  }
  if (is(x, "Rle")) {
    v <- runValue(x)
    if (is.factor(v)) return(enc_vec(as.character(as.vector(x))))
    return(list(class = u("Rle"), values = enc(v)$values, lengths = runLength(x)))
  }
  if (is(x, "DataFrame")) {
    return(list(class = u("DataFrame"), columns = enc_mcols(x), nrow = u(nrow(x))))
  }
  if (is(x, "List")) {
    return(list(class = u("List"), elements = unname(lapply(as.list(x), enc_list_elem))))
  }
  if (is.list(x)) {
    return(list(class = u("List"), elements = unname(lapply(x, enc_list_elem))))
  }
  if (is_atomic_vec(x)) return(enc_vec(x))
  stop("cannot encode value of class ", class(x)[1])
}
