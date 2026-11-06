suppressWarnings(suppressMessages({library(GenomicRanges)}))
args <- commandArgs(TRUE); seed <- as.integer(args[1]); N <- as.integer(args[2])
set.seed(seed)
SI <- Seqinfo(c("b","a","c"), c(40, NA, 30), c(FALSE, NA, FALSE))
gtxt <- function(x) if (length(x)==0) "-" else paste0(as.character(seqnames(x)),"|",start(x),"|",end(x),"|",as.character(strand(x)), collapse=",")
hs <- function(h) if (length(h)==0) "-" else { o <- order(queryHits(h), subjectHits(h)); paste0(queryHits(h)[o]-1L,"-",subjectHits(h)[o]-1L,collapse=",") }
vs <- function(v) if (length(v)==0) "-" else paste(ifelse(is.na(v),"NA",v),collapse=",")
rg <- function(n, maxs=30, maxw=6, minw=0) { if (n==0) return(GRanges(seqinfo=SI)); GRanges(sample(c("b","a","c"),n,TRUE,prob=c(.4,.4,.2)), IRanges(sample(1:maxs,n,TRUE), width=sample(minw:maxw,n,TRUE)), strand=sample(c("+","-","*"),n,TRUE), seqinfo=SI) }
emit <- function(op, par, x, y, res) cat(op, par, gtxt(x), if (is.null(y)) "." else gtxt(y), res, "\n", sep="\t")
try_ <- function(f) tryCatch(f(), error=function(e) "ERR")
for (i in 1:N) {
  n <- sample(0:8,1); m <- sample(0:8,1); x <- rg(n); y <- rg(m); ynz <- rg(m, minw=1); xnz <- rg(n, minw=1)
  for (ig in c(FALSE,TRUE)) {
    P <- as.character(ig)
    r <- reduce(x, ignore.strand=ig, with.revmap=TRUE); emit("reduce", P, x, NULL, paste0(gtxt(r),";",paste(sapply(mcols(r)$revmap,function(v) paste(v-1L,collapse="/")),collapse=",")))
    r <- reduce(x, ignore.strand=ig, min.gapwidth=0, drop.empty.ranges=TRUE); emit("reduce0d", P, x, NULL, gtxt(r))
    emit("range", P, x, NULL, gtxt(range(x, ignore.strand=ig)))
    r <- disjoin(x, ignore.strand=ig, with.revmap=TRUE); emit("disjoin", P, x, NULL, paste0(gtxt(r),";",paste(sapply(mcols(r)$revmap,function(v) paste(v-1L,collapse="/")),collapse=",")))
    emit("gaps", P, x, NULL, try_(function() gtxt(suppressWarnings(gaps(x, ignore.strand=ig)))))
    emit("union", P, x, y, try_(function() gtxt(union(x,y,ignore.strand=ig))))
    emit("intersect", P, x, y, try_(function() gtxt(intersect(x,y,ignore.strand=ig))))
    emit("setdiff", P, x, y, try_(function() gtxt(setdiff(x,y,ignore.strand=ig))))
    for (ty in c("any","start","end","within","equal")) for (mg in c(-1,0,2)) for (mo in c(0,1)) {
      if (ty=="any" && mg>=0 && mo>0) next
      emit("overlaps", paste(P,ty,mg,mo,sep=","), x, y, try_(function() hs(findOverlaps(x,y,type=ty,maxgap=mg,minoverlap=mo,ignore.strand=ig))))
    }
    emit("nearest_ign_arb_in_all", P, xnz, ynz, hs(nearest(xnz,ynz,select="all",ignore.strand=ig)))
    emit("precede", P, xnz, ynz, hs(precede(xnz,ynz,select="all",ignore.strand=ig)))
    emit("follow", P, xnz, ynz, hs(follow(xnz,ynz,select="all",ignore.strand=ig)))
    if (n > 0 && m > 0) { emit("precede_first", P, xnz, ynz, vs(precede(xnz,ynz,ignore.strand=ig)-1L)); emit("follow_first", P, xnz, ynz, vs(follow(xnz,ynz,ignore.strand=ig)-1L)) }
    if (n == m) {
      emit("distance", P, x, y, vs(distance(x,y,ignore.strand=ig)))
      emit("punion", paste(P,"F"), x, y, try_(function() gtxt(punion(x,y,ignore.strand=ig))))
      emit("punion", paste(P,"T"), x, y, try_(function() gtxt(punion(x,y,fill.gap=TRUE,ignore.strand=ig))))
      emit("pintersect", P, x, y, try_(function() gtxt(pintersect(x,y,ignore.strand=ig))))
      emit("psetdiff", P, x, y, try_(function() gtxt(psetdiff(x,y,ignore.strand=ig))))
    }
    if (n > 0) {
      k <- sample(0:8,1); for (fx in c("start","end","center")) emit("resize", paste(P,k,fx,sep=","), x, NULL, gtxt(resize(x,k,fix=fx,ignore.strand=ig)))
      k <- sample(-6:6,1); for (st2 in c(TRUE,FALSE)) for (b in c(FALSE,TRUE)) emit("flank", paste(P,k,st2,b,sep=","), x, NULL, gtxt(flank(x,k,start=st2,both=b,ignore.strand=ig)))
      u <- sample(0:5,1); dd <- sample(0:5,1); emit("promoters", paste(P,u,dd,sep=","), x, NULL, gtxt(promoters(x,u,dd)))
    }
  }
  emit("nearest_arb", "", x, y, vs(nearest(x,y)-1L))
  emit("order", "", x, NULL, vs(order(x)-1L)); emit("orderdec", "", x, NULL, vs(order(x,decreasing=TRUE)-1L))
  emit("rank", "", x, NULL, vs(rank(x,ties.method="min"))); emit("dup", "", x, NULL, vs(as.integer(duplicated(x))))
  emit("unique", "", x, NULL, gtxt(unique(x)))
  emit("trim", "", x, NULL, gtxt(suppressWarnings(trim(x))))
  w <- sample(c(NA,0:30),1); sh <- sample(-2:3,1)
  cv <- if (is.na(w)) coverage(x, shift=sh) else coverage(x, shift=sh, width=w)
  emit("coverage", paste(ifelse(is.na(w),"None",w),sh,sep=","), x, NULL, paste(sapply(names(cv), function(nm) paste0(nm,"=",if(length(cv[[nm]])==0) "-" else paste0(runLength(cv[[nm]]),"x",runValue(cv[[nm]]),collapse="/"))), collapse=";"))
}
