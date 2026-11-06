suppressWarnings(suppressMessages(library(GenomicRanges)))
z <- function(l, v) cat(l, ":", paste(ifelse(is.na(v),"NA",v-1L), collapse=","), "\n")
al <- function(l,h) cat(l, ":", paste(sprintf("%d-%d", queryHits(h)-1L, subjectHits(h)-1L), collapse=" "), "\n")
g <- function(s, st) GRanges("a", IRanges(s, width=2), strand=st)
Q <- g(100, "*")
S <- GRanges("a", IRanges(c(120,80,130,70,60), width=2), strand=c("+","-","*","-","+"))
cat(as.character(S), "\n")
z("precede *", precede(Q,S)); al("precede * all", precede(Q,S,select="all")); z("follow *", follow(Q,S)); al("follow * all", follow(Q,S,select="all"))
S2 <- GRanges("a", IRanges(c(80,70),width=2), strand="-"); z("p only minus left", precede(Q,S2)); z("f only minus left", follow(Q,S2))
S3 <- GRanges("a", IRanges(c(120,130),width=2), strand="-"); z("p only minus right", precede(Q,S3)); z("f only minus right", follow(Q,S3))
S4 <- GRanges("a", IRanges(c(120,70),width=2), strand=c("-","-")); z("p minus both", precede(Q,S4)); z("f minus both", follow(Q,S4)); al("p minus both all", precede(Q,S4,select="all"))
S5 <- GRanges("a", IRanges(c(125,80),width=2), strand=c("-","-")); z("p minus both2", precede(Q,S5)); z("f minus both2", follow(Q,S5))
S6 <- GRanges("a", IRanges(c(125,90),width=2), strand=c("+","-")); z("p +right -left", precede(Q,S6)); z("f +right -left", follow(Q,S6)); al("p all", precede(Q,S6,select="all")); al("f all", follow(Q,S6,select="all"))
S7 <- GRanges("a", IRanges(c(110,70),width=2), strand=c("+","-")); z("p +right(near) -left(far)", precede(Q,S7)); z("f", follow(Q,S7)); 
S8 <- GRanges("a", IRanges(c(140,90),width=2), strand=c("+","-")); z("p +right(far) -left(near)", precede(Q,S8)); z("f", follow(Q,S8)); 
S9 <- GRanges("a", IRanges(c(110,90),width=2), strand=c("*","*")); z("p **", precede(Q,S9)); z("f **", follow(Q,S9));
Qp <- g(100,"+"); z("plus q: p", precede(Qp,S)); z("plus q: f", follow(Qp,S)); al("plus q p all", precede(Qp,S,select="all"))
Qm <- g(100,"-"); z("minus q: p", precede(Qm,S)); z("minus q: f", follow(Qm,S)); al("minus q p all", precede(Qm,S,select="all"));al("minus q f all", follow(Qm,S,select="all"))
cat("nearest * q:\n"); al("nearest *", nearest(Q,S,select="all"))
