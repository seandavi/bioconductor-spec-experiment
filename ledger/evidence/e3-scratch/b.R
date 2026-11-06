source("spec/scratch/h.R")
x <- IRanges(c(1,5,10,20), width=c(3,0,4,1))
p("x",x)
p("shift2", shift(x,2)); p("shiftvec", shift(x,c(1,-1,0,10)))
y<-x[c(1,3)]
p("narrow s", narrow(y,start=2)); 
p("narrow e", narrow(y,end=2))
p("narrow se", narrow(x[3],start=2,end=3))
p("narrow w", narrow(x[3],start=2,width=2))
p("narrow ew", narrow(x[3],end=3,width=2))
p("narrow neg start", narrow(x[3],start=-2))
p("narrow neg end", narrow(x[3],end=-2))
r<-try(narrow(x[1],start=5)); r<-try(narrow(x[2],start=1)); r <- try(narrow(x[1],start=1,end=4))
r<-try(narrow(x[1],start=1,end=3)); p("full",r)
p("narrow zero", try(narrow(x[1],start=3,width=0)))
p("narrow zero2", try(narrow(x[1],start=4,width=0)))
cat("--resize\n")
for (fx in c("start","end","center")) { p(paste("resize 2",fx), resize(x,2,fix=fx)); p(paste("resize 0",fx), resize(x,0,fix=fx));p(paste("resize 5",fx), resize(x,5,fix=fx)) }
p("rs w vec", resize(x, c(1,2,3,4)))
p("rs1 c", resize(IRanges(1,5),2,fix="center")); p("rs1 c", resize(IRanges(1,4),3,fix="center"));p("rs c", resize(IRanges(1,4),2,fix="center")); p("rs c", resize(IRanges(1,3),2,fix="center"));p("rs c", resize(IRanges(2,4),0,fix="center"))
p("rs c", resize(IRanges(1,6),3,fix="center")); p("rs c", resize(IRanges(1,6),1,fix="center"))
p("rs c", resize(IRanges(1,5),4,fix="center")); 
r <- try(resize(x,-1))
cat("--flank\n")
for (b in c(TRUE,FALSE)) for (w in c(2,-2,0)) p(paste("flank",b,w), flank(x,w,start=b))
p("flank both", flank(x,2,both=TRUE)); p("flank both end", flank(x,2,start=FALSE,both=TRUE));p("flank both neg", flank(x,-2,both=TRUE));
p("flank both 0", flank(x,0,both=TRUE));p("flank both 1", flank(x,1,both=TRUE))
cat("--promoters\n")
p("prom", promoters(x)); p("prom 3 1", promoters(x, upstream=3, downstream=1)); p("prom 0 0", promoters(x,0,0));p("prom 0 5", promoters(x,0,5))
p("prom 2 0", promoters(x,2,0)); 
