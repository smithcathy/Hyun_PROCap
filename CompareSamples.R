suppressMessages(library(Rcpp))
suppressMessages(library(dplyr))

tic=proc.time()
args=as.numeric(commandArgs(trailingOnly = T))
#args[1] is distance, args[2] is id

data=read.table(file("stdin"), stringsAsFactors = F)
data=select(data,V1,V2,V5,V7)
data=filter(data, V1 %in% c(1:22,"X"))
#data=filter(data, V1 %in% c(1:22,"X","Y"))
data=data %>% mutate(V1=replace(V1, V1=="X", "23"))
#data=data %>% mutate(V1=replace(V1, V1=="Y", "24"))
data$V1=as.numeric(data$V1)
maxLengthByChr=aggregate(data$V2,by=list(data$V1),max)$x
totalLength=sum(as.numeric(maxLengthByChr))

sourceCpp("CompareSamples.cpp")
r=GenomeWideCompare(data$V1,data$V2, data$V5, data$V7,args[1], totalLength)
reads=as.data.frame(r$reads)
reads=filter(reads, V2!=0)
file=paste("SamplesToCompare/ReadsAtLoci",args[2],args[1],sep="_")
write.table(reads,file,row.names = F, col.names = F)
toc=proc.time()
print(toc-tic)