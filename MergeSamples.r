suppressMessages(library(Rcpp))
suppressMessages(library(dplyr))
suppressMessages(library(hexbin))

#tic=proc.time()
args=commandArgs(trailingOnly = T)
#args[1] is filename1, args[2] is filename2
id1=strsplit(args[1],"_")[[1]][2]
id2=strsplit(args[2],"_")[[1]][2]
dist=strsplit(args[1],"_")[[1]][3]

data1=as.matrix(read.table(args[1],stringsAsFactors = F, header = T, colClasses = c("numeric","numeric","numeric")))
data2=as.matrix(read.table(args[2],stringsAsFactors = F, header = T, colClasses = c("numeric","numeric","numeric")))

sourceCpp("MergeSamples.cpp")
r=MergeSamples(data1,data2)
reads=as.data.frame(r$mergedData)
reads=filter(reads, V1!=0)
outfile=paste("SamplesToCompare/ScatterPlots",paste("GW",id1,"V",id2,"_",dist,"Truncated",".pdf",sep = ""),sep="/")
title=paste("Counts:",3000000000-r$counts[1],r$counts[2],"\n",r$counts[3],r$counts[4],sep = " ")
pdf(outfile)
#x=gplot.hexbin(hexbin(reads$V3,reads$V4),xlab=id1,ylab = id2,main=title,colorcut = 100000,
               #colramp = function(n) LinGray(n,beg = 90,end = 10))
x=gplot.hexbin(hexbin(reads$V3,reads$V4),xlab=id1,ylab = id2,main=title,colorcut = 10000,maxcnt = 2000,mincnt = 10,
               colramp = function(n) BTC(n,beg = 30,end = 240))
dev.off()
#toc=proc.time()
#print(toc-tic)