source("spec/scratch/h.R")
z <- function(l, v) cat(l, ":", paste(ifelse(is.na(v),"NA",v-1L), collapse=","), "\n")
al <- function(l,h) cat(l, ":", paste(sprintf("%d-%d", queryHits(h)-1L, subjectHits(h)-1L), collapse=" "), "\n")
Q <- IRanges(6,7)
S <- IRanges(c(1,10),c(3,12)); al("tie", nearest(Q,S,select="all")); z("arb", nearest(Q,S))
S <- IRanges(c(10,1),c(12,3)); al("tie sw", nearest(Q,S,select="all")); z("arb", nearest(Q,S))
S <- IRanges(c(1,10,10,1),c(3,12,13,2)); al("tie 4", nearest(Q,S,select="all")); z("arb", nearest(Q,S))
S <- IRanges(c(1,9),c(3,12)); al("nontie", nearest(Q,S,select="all"))
S <- IRanges(c(1,8),c(3,12)); al("nontie", nearest(Q,S,select="all"))
# random check of nearest all vs min-distance
set.seed(1)
bad <- 0
for (i in 1:2000) {
  n <- sample(1:6,1); m <- sample(1:6,1)
  ss <- sample(1:30,n,TRUE); sw <- sample(0:5,n,TRUE); qs <- sample(1:30,m,TRUE); qw <- sample(0:4,m,TRUE)
  S <- IRanges(ss,width=sw); Q <- IRanges(qs,width=qw)
  h <- nearest(Q,S,select="all")
  d <- outer(seq_along(Q), seq_along(S), Vectorize(function(a,b) max(0, max(start(Q)[a]-end(S)[b], start(S)[b]-end(Q)[a]) - 1)))
  exp <- which(d == apply(d,1,min)[row(d)], arr.ind=TRUE)
  exp <- exp[order(exp[,1],exp[,2]),,drop=FALSE]
  got <- cbind(queryHits(h), subjectHits(h)); got <- got[order(got[,1],got[,2]),,drop=FALSE]
  if (!identical(unname(exp), unname(got))) { bad <- bad+1; if (bad<=4) { cat("S:", s(S), " Q:", s(Q), "\n"); print(unname(exp)); print(got) } }
}
cat("bad", bad, "\n")
