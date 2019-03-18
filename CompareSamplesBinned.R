library(Rcpp)
library(dplyr)
library(feather)

tic=proc.time()
args=as.numeric(commandArgs(trailingOnly = T))
#args[1] is distance, args[2] is binSize, args[3] is id

data=read.table(file("stdin"), stringsAsFactors = F)
data=select(data,V1,V2,V5,V7)
data=filter(data, V1 %in% c(1:22,"X","Y"))
data=data %>% mutate(V1=replace(V1, V1=="X", "23"))
data=data %>% mutate(V1=replace(V1, V1=="Y", "24"))
data$V1=as.numeric(data$V1)
str(data)
chrUsed=unique(data$V1)
maxLengthByChr=aggregate(data$V2,by=list(data$V1),max)$x
print(maxLengthByChr)
totalLength=sum(as.numeric(maxLengthByChr))
maxChrLength=max(data$V2)

sourceCpp("CompareSamplesBinned.cpp")
r=GenomeWideDistBin(data$V1,data$V2, data$V5, data$V7,floor(args[1]/args[2]),maxChrLength,chrUsed,maxLengthByChr,
                    args[2], ceiling(totalLength/args[2]))
str(r)
file=paste("ReadsAtLoci",args[3],sep="_")
write_feather(as.data.frame(r$reads),file)
toc=proc.time()
print(toc-tic)