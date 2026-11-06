source("spec/scratch/h.R")
q <- function(l, e) { r <- try(e, silent=TRUE); if (inherits(r,"try-error")) cat(l,": ERROR\n") else p(l, r) }
x <- IRanges(10,13)
q("end0", narrow(x,end=0)); q("end -4", narrow(x,end=-4)); q("end -5", narrow(x,end=-5)); q("start -5", narrow(x,start=-5)); q("start -4", narrow(x,start=-4))
q("start 5", narrow(x,start=5)); q("start 5 w0", narrow(x,start=5,width=0)); q("start 4 end 3", narrow(x,start=4,end=3)); q("start 4 end 2", narrow(x,start=4,end=2))
q("end 1 w 0", narrow(x,end=1,width=0)); q("end 0 w 0", narrow(x,end=0,width=0)); q("start1 w5", narrow(x,start=1,width=5)); q("end 4 w 5", narrow(x,end=4,width=5))
q("start end w", narrow(x,start=1,end=3,width=3))
q("zero start1", narrow(IRanges(5,4),start=1)); q("zero end0", narrow(IRanges(5,4),end=0)); q("zero start 2", narrow(IRanges(5,4),start=2))
q("start -1", narrow(x,start=-1)); q("end -1 start -3", narrow(x,start=-3,end=-2))
