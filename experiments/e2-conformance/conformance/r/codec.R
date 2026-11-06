# R codec for the conformance format (see ../FORMAT.md).
#   encode(x)          R object -> tree (nested lists; scalars are length-1 atomics)
#   decode(tree, fx)   tree -> R object; fx = named list of decoded fixtures
#   to_json(tree)      tree -> JSON text (stable key order, arrays of scalars inline)
#   from_json(text)    JSON text -> tree
#   tree_equal(a, b)   deep equality with numeric tolerance

suppressPackageStartupMessages(suppressWarnings({
  library(GenomicRanges)
  library(jsonlite)
}))

# ---- encode ---------------------------------------------------------------

enc_items <- function(x, tag) {
  items <- as.list(unname(x))
  if (tag == "$dbl") {
    items <- lapply(items, function(v) {
      if (is.nan(v)) "NaN" else if (is.na(v)) NA_real_
      else if (v == Inf) "Inf" else if (v == -Inf) "-Inf" else v
    })
  }
  items
}

enc_atomic <- function(x) {
  tag <- if (is.integer(x)) "$int" else if (is.double(x)) "$dbl"
         else if (is.character(x)) "$chr" else if (is.logical(x)) "$lgl"
         else stop("cannot encode atomic of class ", class(x)[1])
  out <- setNames(list(enc_items(x, tag)), tag)
  if (!is.null(names(x))) out$names <- enc_items(names(x), "$chr")
  out
}
enc_int <- function(x) enc_atomic(as.integer(x))
enc_chr <- function(x) enc_atomic(as.character(x))
enc_names <- function(x) if (is.null(names(x))) NULL else enc_chr(names(x))

enc_dframe <- function(df, n) {
  cols <- if (is.null(df)) list() else
    setNames(lapply(seq_len(ncol(df)), function(j) encode(df[[j]])), colnames(df))
  if (length(cols) == 0) cols <- setNames(list(), character(0))
  list(`$type` = "DFrame", nrow = n, columns = cols)
}

enc_seqinfo <- function(si) {
  list(`$type` = "Seqinfo", seqnames = enc_chr(seqnames(si)),
       seqlengths = enc_int(seqlengths(si)),
       isCircular = enc_atomic(unname(isCircular(si))),
       genome = enc_chr(unname(genome(si))))
}

encode <- function(x) {
  if (is.null(x)) return(NULL)
  if (is(x, "GRanges")) {
    list(`$type` = "GRanges",
         seqnames = enc_chr(as.character(seqnames(x))),
         start = enc_int(start(x)), width = enc_int(width(x)),
         strand = enc_chr(as.character(strand(x))),
         names = enc_names(x),
         mcols = enc_dframe(mcols(x, use.names = FALSE), length(x)),
         seqinfo = enc_seqinfo(seqinfo(x)))
  } else if (is(x, "IntegerRanges")) {
    list(`$type` = "IRanges",
         start = enc_int(start(x)), width = enc_int(width(x)),
         names = enc_names(x),
         mcols = enc_dframe(mcols(x, use.names = FALSE), length(x)))
  } else if (is(x, "Hits")) {
    list(`$type` = "Hits", from = enc_int(queryHits(x)), to = enc_int(subjectHits(x)),
         nLnode = nLnode(x), nRnode = nRnode(x))
  } else if (is(x, "Seqinfo")) {
    enc_seqinfo(x)
  } else if (is(x, "RleList")) {
    els <- setNames(lapply(as.list(x), encode), names(x))
    list(`$type` = "RleList", elements = els)
  } else if (is(x, "Rle")) {
    list(`$type` = "Rle", values = encode(runValue(x)), lengths = enc_int(runLength(x)))
  } else if (is(x, "IntegerList")) {
    list(`$type` = "IntegerList",
         elements = lapply(as.list(x), function(v) as.list(as.integer(v))),
         names = if (is.null(names(x))) NULL else enc_chr(names(x)))
  } else if (is(x, "DataFrame")) {
    enc_dframe(x, nrow(x))
  } else if (is.factor(x)) {
    list(`$type` = "Factor", levels = enc_chr(levels(x)), codes = enc_int(as.integer(x)))
  } else if (is.atomic(x)) {
    enc_atomic(x)
  } else stop("cannot encode object of class ", class(x)[1])
}

# ---- JSON text ------------------------------------------------------------

TAGS <- c("$int", "$dbl", "$chr", "$lgl")

json_scalar <- function(v) {
  if (is.null(v) || (length(v) == 1 && is.na(v) && !is.nan(v))) return("null")
  if (is.character(v)) return(as.character(jsonlite::toJSON(v, auto_unbox = TRUE)))
  if (is.logical(v)) return(if (v) "true" else "false")
  if (is.integer(v)) return(as.character(v))
  if (is.double(v)) {
    if (is.nan(v)) return("\"NaN\"")
    s <- sprintf("%.15g", v)
    return(if (grepl("^-?[0-9]+$", s)) paste0(s, ".0") else s)
  }
  stop("bad scalar")
}

is_scalar <- function(v) is.null(v) || (is.atomic(v) && length(v) == 1)

to_json <- function(tree, indent = 0) {
  pad <- strrep("  ", indent)
  pad1 <- strrep("  ", indent + 1)
  if (is_scalar(tree)) return(json_scalar(tree))
  if (!is.list(tree)) stop("bad tree node of class ", class(tree)[1])
  if (length(tree) == 0) return(if (is.null(names(tree))) "[]" else "{}")
  if (is.null(names(tree))) {
    if (all(vapply(tree, is_scalar, TRUE)))
      return(paste0("[", paste(vapply(tree, json_scalar, ""), collapse = ", "), "]"))
    items <- vapply(tree, to_json, "", indent = indent + 1)
    return(paste0("[\n", paste0(pad1, items, collapse = ",\n"), "\n", pad, "]"))
  }
  if (any(names(tree) %in% TAGS)) {  # tagged vector: one line
    return(paste0("{", paste0(vapply(names(tree), json_scalar, ""), ": ",
                              vapply(tree, to_json, ""), collapse = ", "), "}"))
  }
  keys <- vapply(names(tree), json_scalar, "")
  vals <- vapply(seq_along(tree), function(i) to_json(tree[i][[1]], indent + 1), "")
  paste0("{\n", paste0(pad1, keys, ": ", vals, collapse = ",\n"), "\n", pad, "}")
}

from_json <- function(text) jsonlite::fromJSON(text, simplifyVector = FALSE)
read_json <- function(path) from_json(paste(readLines(path, warn = FALSE), collapse = "\n"))

# ---- equality -------------------------------------------------------------

tree_equal <- function(a, b) {
  if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
  if (is.list(a) || is.list(b)) {
    if (!(is.list(a) && is.list(b))) return(FALSE)
    if (length(a) != length(b)) return(FALSE)
    if (length(a) == 0) return(TRUE)
    if (!identical(names(a), names(b))) return(FALSE)
    return(all(vapply(seq_along(a), function(i) tree_equal(a[i][[1]], b[i][[1]]), TRUE)))
  }
  if (is.numeric(a) && is.numeric(b))
    return(isTRUE(abs(a - b) <= 1e-9 * max(1, abs(b))))
  identical(a, b)
}

# ---- decode ---------------------------------------------------------------

dec_items <- function(items, tag) {
  proto <- switch(tag, `$int` = integer(1), `$dbl` = numeric(1),
                  `$chr` = character(1), `$lgl` = logical(1))
  na <- switch(tag, `$int` = NA_integer_, `$dbl` = NA_real_,
               `$chr` = NA_character_, `$lgl` = NA)
  vapply(items, function(v) {
    if (is.null(v)) na else if (tag == "$dbl" && is.character(v)) as.numeric(v)
    else if (tag == "$int") as.integer(v) else v
  }, proto, USE.NAMES = FALSE)
}

dec_shorthand <- function(v) {
  items <- if (is.list(v)) v else list(v)
  if (length(items) == 0) stop("ambiguous empty array: use a tag")
  nn <- Filter(Negate(is.null), items)
  if (length(nn) == 0) return(rep(NA, length(items)))
  kinds <- unique(vapply(nn, function(e) class(e)[1], ""))
  if (any(vapply(nn, is.list, TRUE))) stop("bad shorthand value")
  if (length(kinds) == 1 && kinds == "character") return(dec_items(items, "$chr"))
  if (length(kinds) == 1 && kinds == "logical") return(dec_items(items, "$lgl"))
  if (all(kinds %in% c("integer", "numeric")))
    return(dec_items(items, if ("numeric" %in% kinds) "$dbl" else "$int"))
  stop("mixed kinds in shorthand array")
}

dec_vec <- function(tree) {
  tag <- intersect(names(tree), TAGS)
  x <- dec_items(tree[[tag]], tag)
  if (!is.null(tree$names)) names(x) <- dec_items(tree$names, "$chr")
  x
}

dec_ranges <- function(tree, fx) {
  g <- function(k) if (is.null(tree[[k]])) NULL else decode(tree[[k]], fx)
  x <- IRanges(start = g("start"), end = g("end"), width = g("width"), names = g("names"))
  if (!is.null(tree$mcols)) mcols(x) <- decode(tree$mcols, fx)
  x
}

decode <- function(tree, fx = list()) {
  if (is.null(tree)) return(NULL)
  if (!is.list(tree)) return(dec_shorthand(tree))
  if (is.null(names(tree))) return(dec_shorthand(tree))
  if (!is.null(tree$`$ref`)) {
    if (is.null(fx[[tree$`$ref`]])) stop("unknown fixture ", tree$`$ref`)
    return(fx[[tree$`$ref`]])
  }
  if (any(names(tree) %in% TAGS)) return(dec_vec(tree))
  type <- tree$`$type`
  if (is.null(type)) stop("object without $type")
  g <- function(k) if (is.null(tree[[k]])) NULL else decode(tree[[k]], fx)
  suppressWarnings(switch(type,
    IRanges = dec_ranges(tree, fx),
    Seqinfo = Seqinfo(seqnames = g("seqnames"), seqlengths = g("seqlengths"),
                      isCircular = g("isCircular"), genome = g("genome")),
    GRanges = {
      sn <- g("seqnames")
      si <- if (!is.null(tree$seqinfo)) g("seqinfo") else Seqinfo(unique(sn))
      gr <- GRanges(seqnames = factor(sn, levels = seqlevels(si)),
                    ranges = dec_ranges(tree, fx),
                    strand = if (is.null(tree$strand)) "*" else g("strand"),
                    seqinfo = si)
      gr
    },
    DFrame = {
      cols <- lapply(tree$columns, decode, fx = fx)
      n <- if (!is.null(tree$nrow)) tree$nrow else if (length(cols)) length(cols[[1]]) else 0L
      df <- S4Vectors::make_zero_col_DFrame(n)
      for (nm in names(cols)) df[[nm]] <- cols[[nm]]
      df
    },
    Hits = Hits(from = g("from"), to = g("to"), nLnode = g("nLnode"),
                nRnode = g("nRnode"), sort.by.query = FALSE),
    Rle = Rle(g("values"), g("lengths")),
    RleList = RleList(lapply(tree$elements, decode, fx = fx)),
    IntegerList = {
      l <- lapply(tree$elements, function(v) as.integer(unlist(v)))
      if (!is.null(tree$names)) names(l) <- g("names")
      IntegerList(l)
    },
    Factor = {
      lv <- g("levels"); cd <- g("codes")
      factor(lv[cd], levels = lv)
    },
    stop("unknown $type ", type)))
}
