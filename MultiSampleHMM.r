suppressMessages(library(Rcpp))

tic=proc.time()
args=commandArgs(trailingOnly = T)
print(args)
binSize=as.numeric(args[1])
windowSize=as.numeric(args[2])
dist=as.numeric(args[3])
files=args[4:length(args)]

sourceCpp("MultiSampleHMM.cpp")
r=MultiSampleHMM(binSize,windowSize,dist,files)
print(r[1])
r=r[-1]
#print(r)
#str(r)
attributes(r)=list(row.names=c(NA_integer_,5), class="data.frame", names=c(1:length(r)))
str(r)
r=t(r)
str(r)
outfile="MergedSamples"
write.table(r,outfile,row.names = F, col.names = F)
toc=proc.time()
print(toc-tic)