library(Rcpp)

chrLength=c(248956422,242193529,198295559,19021455,181538259,170805979,159345973,145138636,138394717,133797422,135086622,
            133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,
            156040895,57227415) #stores length of each chromosome (from wiki)

extractLoci=function(loci){
  if (length(strsplit(loci,":")[[1]])==1){
    return(c(0,chrLength[as.numeric(strsplit(loci,":")[[1]][1])])) #returns length of full chromosome
  }else{
    startAndEnd=strsplit(loci,":")[[1]][2] #splits at semi-colon to separate chromosome and loci
    start=as.numeric(strsplit(startAndEnd,"-")[[1]][1]) #takes start site
    end=as.numeric(strsplit(startAndEnd,"-")[[1]][2]) #takes end site
    return(c(start,end-start))
  }
}

setwd("Documents/HyunLab")
seq=read.table("48646_R1_cap.slfn5.bed",stringsAsFactors = F)
loci="17:33550000-33590000"
binSize=10
interval=400
trans=matrix(c(.999,.001,.99,.01),nrow = 2, ncol = 2) # (i,j) means it came from state j and is moving to state i
sourceCpp("PROcap_BestDist.cpp")
tic=proc.time()
seqResults=PROcapBestDist(seq$V1, seq$V2, seq$V5, seq$V7, extractLoci(loci)[1], extractLoci(loci)[2], binSize, 
                          interval,log(trans))
toc=proc.time()
toc-tic
sourceCpp("PROcap_scaled.cpp")
tic=proc.time()
seqResults=PROcapScaled(seq$V1, seq$V2, seq$V5, seq$V7, extractLoci(loci)[1], extractLoci(loci)[2], binSize, windowSize,dist,log(trans))
toc=proc.time()
toc-tic

fakeseq=seq
fakeseq$V5=c(rpois(492,1),50,rpois(143,1))
fakeseq$V7=c(rpois(422,1),50,rpois(213,1))
#expect change of state at i=250
sourceCpp("PROCap_scaled.cpp")
tic=proc.time()
fakeseqResults=PROcapScaled(fakeseq$V1, fakeseq$V2, fakeseq$V5, fakeseq$V7, extractLoci(loci)[1], extractLoci(loci)[2], 
                            binSize, windowSize,dist,log(trans))
toc=proc.time()
toc-tic

chr17=read.table("48646_R1_cap.chr17.bed",stringsAsFactors = F)
loci="17:"
sourceCpp("PROcap_BestDist.cpp")
tic=proc.time()
chr17Results=PROcapBestDist(chr17$V1, chr17$V2, chr17$V5, chr17$V7, extractLoci(loci)[1], extractLoci(loci)[2], binSize, 
                          interval,log(trans))
toc=proc.time()
toc-tic
sourceCpp("PROCap_scaled.cpp")
tic=proc.time()
chr17Results=PROcapScaled(chr17$V1, chr17$V2, chr17$V5, chr17$V7, extractLoci(loci)[1], extractLoci(loci)[2], binSize, 
                    windowSize,dist,log(trans))
toc=proc.time()
toc-tic


getFile = function(dir,loci){
  file=strsplit(dir,"/")[[1]][str_count(dir,"/")+1]
  path=substr(dir,1,nchar(dir)-nchar(file))
  command=pipe(paste(paste("cd",path,sep = " "),
                     paste("tabix -pbed",file,sep = " "),sep = " | "))
  result=capture.output(command)
  return(result)
}



findTFBS = function(dir, loci, dist, lambda, binSize){
  seq=getFile(dir,loci) #formats file in shell and then imports
  sourceCpp("PROcap.cpp")
  TFBS=PROcap(select(seq,V1),select(seq,V2),select(seq,V5),loci,dist,lamda,binSize)
}


get_file(paste(getwd(),"48651_R1_cap.clip.ACAGTG.hs37d5.bowtie.bed.gz",sep = "/"),"abc")

#tabix input bed file chr17:
#rscript [script.r]
#/dex/stdin open file fopen("dev.stdin") or fstream or ifstream recc ifstream

babyseq=seq[1:13,]
babyseq$V5=c(rpois(3,5),30,rpois(9,5))
babyseq$V7=c(30,rpois(12,5))
#expect change at i=6
loci="17:33568000-33568150"
binSize=10
sourceCpp("PROCap_scaled.cpp")
tic=proc.time()
babyseqResults=PROcapScaled(babyseq$V1, babyseq$V2, babyseq$V5, babyseq$V7, extractLoci(loci)[1], extractLoci(loci)[2], 
                            binSize, windowSize,dist,log(trans))
toc=proc.time()
toc-tic
