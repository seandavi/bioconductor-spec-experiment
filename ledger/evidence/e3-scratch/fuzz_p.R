suppressWarnings(suppressMessages({library(IRanges)}))
args <- commandArgs(TRUE); seed <- as.integer(args[1]); N <- as.integer(args[2])
set.seed(seed)
rs <- function(x) if (length(x)==0) "-" else paste0(start(x),":",end(x),collapse=",")
hs <- function(h) if (length(h)==0) "-" else { o <- order(queryHits(h), subjectHits(h)); paste0(queryHits(h)[o]-1L,"-",subjectHits(h)[o]-1L,collapse=",") }
vs <- function(v) if (length(v)==0) "-" else paste(ifelse(is.na(v),"NA",v),collapse=",")
rr <- function(n, maxs=25, maxw=6, zero=TRUE) IRanges(start=sample(1:maxs,n,TRUE), width=sample(if(zero) 0:maxw else 1:maxw,n,TRUE))
emit <- function(op, par, x, y, res) cat(op, par, rs(x), if (is.null(y)) "." else rs(y), res, sep="\t", fill=FALSE, "\n")
try_ <- function(f) tryCatch(f(), error=function(e) "ERR")
for (i in 1:N) {
  n <- sample(1:5,1); m <- n; x <- rr(n, 14, 5); y <- rr(m, 14, 5)
  for (mg in c(0,1,2,3)) for (de in c(FALSE,TRUE)) {
    r <- reduce(x, min.gapwidth=mg, drop.empty.ranges=de, with.revmap=TRUE)
    emit("reduce", paste0(mg,",",de), x, NULL, paste0(rs(r),";",paste(sapply(mcols(r)$revmap,function(v) paste(v-1L,collapse="/")),collapse=",")))
  }
  emit("range", "", x, NULL, if (n) rs(range(x)) else "-")
  st <- sample(c(NA,-2:30),1); en <- sample(c(NA,-2:30),1)
  r <- try_(function() { a <- list(x); if(!is.na(st)) a$start <- st; if(!is.na(en)) a$end <- en; rs(do.call(gaps, a)) })
  emit("gaps", paste0(ifelse(is.na(st),"None",st),",",ifelse(is.na(en),"None",en)), x, NULL, r)
  d <- disjoin(x, with.revmap=TRUE); emit("disjoin", "", x, NULL, paste0(rs(d),";",paste(sapply(mcols(d)$revmap,function(v) paste(v-1L,collapse="/")),collapse=",")))
  emit("union", "", x, y, rs(union(x,y))); emit("intersect", "", x, y, rs(intersect(x,y))); emit("setdiff", "", x, y, rs(setdiff(x,y)))
  if (n == m) {
    emit("punion", "FALSE", x, y, try_(function() rs(punion(x,y))))
    emit("punion", "TRUE", x, y, try_(function() rs(punion(x,y,fill.gap=TRUE))))
    for (re in c("none","max.start","start.x")) emit("pintersect", re, x, y, try_(function() rs(pintersect(x,y,resolve.empty=re))))
    emit("psetdiff", "", x, y, try_(function() rs(psetdiff(x,y))))
    emit("distance", "", x, y, vs(distance(x,y)))
  }
  w <- sample(c(NA,0:30),1); sh <- sample(-3:3,1)
  r <- try_(function() { cv <- if (is.na(w)) coverage(x, shift=sh) else coverage(x, shift=sh, width=w); if (length(cv)==0) "-" else paste0(runLength(cv),"x",runValue(cv),collapse=",") })
  emit("coverage", paste0(ifelse(is.na(w),"None",w),",",sh), x, NULL, r)
  for (ty in c("any","start","end","within","equal")) for (mg in c(-1,0,1,3)) for (mo in c(0,1,2)) {
    if (ty=="any" && mg>=0 && mo>0) next
    emit("overlaps", paste0(ty,",",mg,",",mo), x, y, try_(function() hs(findOverlaps(x,y,type=ty,maxgap=mg,minoverlap=mo))))
  }
  emit("nearest", "", x, y, hs(nearest(x,y,select="all")))
  emit("nearest_self", "", x, NULL, if (n>0) hs(nearest(x,select="all")) else "-")
  emit("precede", "", x, y, hs(precede(x,y,select="all"))); emit("follow", "", x, y, hs(follow(x,y,select="all")))
  emit("precede_first", "", x, y, vs(precede(x,y)-1L)); emit("follow_first", "", x, y, vs(follow(x,y)-1L))
  emit("order", "", x, NULL, vs(order(x)-1L)); emit("orderdec", "", x, NULL, vs(order(x,decreasing=TRUE)-1L))
  for (tm in c("average","min","max","first")) emit("rank", tm, x, NULL, vs(rank(x,ties.method=tm)))
  emit("dup", "", x, NULL, vs(as.integer(duplicated(x)))); emit("dupl", "", x, NULL, vs(as.integer(duplicated(x,fromLast=TRUE))))
  if (n > 0) {
    k <- sample(-5:8,1); emit("shift", k, x, NULL, rs(shift(x,k)))
    k <- sample(0:8,1); for (fx in c("start","end","center")) emit("resize", paste0(k,",",fx), x, NULL, rs(resize(x,k,fix=fx)))
    k <- sample(-6:6,1); for (st2 in c(TRUE,FALSE)) for (b in c(FALSE,TRUE)) emit("flank", paste0(k,",",st2,",",b), x, NULL, rs(flank(x,k,start=st2,both=b)))
    u <- sample(0:5,1); dd <- sample(0:5,1); emit("promoters", paste0(u,",",dd), x, NULL, rs(promoters(x,u,dd)))
    a <- sample(c(NA,-7:7),1); b <- sample(c(NA,-7:7),1); wd <- sample(c(NA,0:6),1)
    r <- try_(function() { ar <- list(x); if(!is.na(a)) ar$start <- a; if(!is.na(b)) ar$end <- b; if(!is.na(wd)) ar$width <- wd; rs(do.call(narrow, ar)) })
    emit("narrow", paste0(ifelse(is.na(a),"None",a),",",ifelse(is.na(b),"None",b),",",ifelse(is.na(wd),"None",wd)), x, NULL, r)
  }
}
