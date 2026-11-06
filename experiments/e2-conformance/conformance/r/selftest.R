# Codec roundtrip on hand-built objects: encode(decode(encode(x))) == encode(x), via JSON text.
args <- commandArgs(FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
source(file.path(here, "codec.R"))

si <- Seqinfo(c("chr2", "chr1"), c(100L, NA), c(FALSE, TRUE), c("hg", NA))
gr <- GRanges(factor(c("chr1", "chr2", "chr1"), levels = seqlevels(si)),
              IRanges(c(1, 5, 9), width = c(3, 0, 2), names = c("a", "b", NA)),
              strand = c("+", "-", "*"), seqinfo = si, score = c(1.5, NA, Inf))
ir <- IRanges(c(1, 5), width = c(3, 0))
mcols(ir)$x <- c("u", NA)
objs <- list(
  ir = ir, gr = gr, si = si,
  rv = reduce(IRanges(c(1, 2, 10), width = 3), with.revmap = TRUE),
  hits = findOverlaps(ir, ir), cov = coverage(ir),
  grcov = coverage(gr), int = c(a = 1L, b = NA), dbl = c(1.5, NaN, -Inf, NA),
  chr = c("x", NA), lgl = c(TRUE, NA), fac = factor(c("a", "b", NA)),
  seqn = seqnames(gr), nul = NULL, empty = integer(0),
  dfl = DataFrame(groups = IntegerList(1:2, integer(0), 3L), f = factor(c("lo", NA, "hi"), levels = c("lo", "hi", "unused"))))
bad <- 0
for (nm in names(objs)) {
  t1 <- from_json(to_json(encode(objs[[nm]])))
  t2 <- from_json(to_json(encode(decode(t1))))
  ok <- tree_equal(t1, t2)
  cat(sprintf("%-6s %s\n", nm, if (ok) "ok" else "MISMATCH"))
  if (!ok) bad <- bad + 1
}
quit(status = as.integer(bad > 0))
