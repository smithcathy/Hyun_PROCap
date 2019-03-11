library(Rcpp)
library(dplyr)

tic=proc.time()
args=as.numeric(commandArgs(trailingOnly = T))
#args[1] is interval

data=read.table(file("stdin"), stringsAsFactors = F)
data=filter(data, V1 %in% c(1:22,"X","Y"))
data=data %>% mutate(V1=replace(V1, V1=="X", "23"))
data=data %>% mutate(V1=replace(V1, V1=="Y", "24"))
data$V1=as.numeric(data$V1)
str(data)
chrUsed=unique(data$V1)
maxLengthByChr=aggregate(data$V2,by=list(data$V1),max)$x
maxChrLength=max(data$V2)

sourceCpp("GenomeWideTest.cpp")
print("cpp file loaded!")
r=GenomeWideDist(data$V1,data$V2, data$V5, data$V7,args[1],maxChrLength,chrUsed,maxLengthByChr)
print("exited cpp")
print(r$maxDist)
print(r$forLambda)
print(r$revLambda)

pdf("GWReadsPerDistLog.pdf")
plot(r$readsPerDistLog)
dev.off()

pdf("GWReadsPerDistSqrt.pdf")
plot(r$readsPerDistSqrt)
dev.off()
toc=proc.time()
print(toc-tic)