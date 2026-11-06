suppressWarnings(source("spec/scratch/h.R"))
set.seed(2)
bad <- 0
for (i in 1:3000) {
  n <- sample(1:6,1); m <- sample(1:3,1)
  ss <- sample(1:30,n,TRUE); sw <- sample(1:5,n,TRUE); qs <- sample(1:30,m,TRUE); qw <- sample(1:4,m,TRUE)
  stq <- sample(c("+","-","*"),1)
  S <- GRanges("a", IRanges(ss,width=sw), strand="+"); Q <- GRanges("a", IRanges(qs,width=qw), strand=stq)
  strand(S) <- sample(c("+","-","*"), n, TRUE)
  h <- nearest(Q,S,select="all")
  got <- cbind(queryHits(h), subjectHits(h)); got <- got[order(got[,1],got[,2]),,drop=FALSE]
  expl <- list()
  for (a in seq_along(Q)) {
    ok <- which(as.character(strand(S)) == as.character(strand(Q))[a] | as.character(strand(S))=="*" | as.character(strand(Q))[a]=="*")
    if (!length(ok)) next
    d <- sapply(ok, function(b) max(0, max(start(Q)[a]-end(S)[b], start(S)[b]-end(Q)[a]) - 1))
    sel <- ok[d == min(d)]
    if (min(d) > 0) {
      up <- if (as.character(strand(Q))[a] == "-") end(S)[sel] < start(Q)[a] else start(S)[sel] > end(Q)[a]
      if (FALSE) sel <- sel[up] 
    }
    expl[[length(expl)+1]] <- cbind(a, sel)
  }
  exp <- if (length(expl)) do.call(rbind, expl) else matrix(integer(),0,2)
  exp <- exp[order(exp[,1],exp[,2]),,drop=FALSE]
  if (!identical(unname(exp), unname(got))) { bad <- bad+1; if (bad<=4) { cat("S:", as.character(S), " Q:", as.character(Q), "\n"); print(unname(exp)); print(got) } }
}
cat("bad", bad, "\n")
