library(Rcpp)
library(VennDiagram)

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

setwd("/Users/cathylicioussmith/Documents/HyunLab")
binSize=c(10,20)
dists=seq(100,400,by=10)
windowSize=500
trans=matrix(c(.999,.001,.9,.1),nrow = 2, ncol = 2) # (i,j) means it came from state j and is moving to state i

seqSLFN_48646=read.table("48646_R1_cap.slfn5.bed",stringsAsFactors = F)
seq17_48646=read.table("48646_R1_cap.chr17.bed",stringsAsFactors = F)
seqSLFN_51644=read.table("51644_R1_cap.slfn5.bed",stringsAsFactors = F)
seq17_51644=read.table("51644_R1_cap.chr17.bed",stringsAsFactors = F)
seqSLFN_48641=read.table("48641_R1_cap.slfn5.bed",stringsAsFactors = F)
seq17_48641=read.table("48641_R1_cap.chr17.bed",stringsAsFactors = F)
seqSLFN_48651=read.table("48651_R1_cap.slfn5.bed",stringsAsFactors = F)
seq17_48651=read.table("48651_R1_cap.chr17.bed",stringsAsFactors = F)

loci="17:33550000-33590000"
for (bin in binSize){
  for (dist in dists){
    slfnResults_48646=PROcapScaled(seqSLFN_48646$V1, seqSLFN_48646$V2, seqSLFN_48646$V5, seqSLFN_48646$V7, extractLoci(loci)[1], extractLoci(loci)[2], bin, windowSize,dist,log(trans))
    slfnResults_51644=PROcapScaled(seqSLFN_51644$V1, seqSLFN_51644$V2, seqSLFN_51644$V5, seqSLFN_51644$V7, extractLoci(loci)[1], extractLoci(loci)[2], bin, windowSize,dist,log(trans))
    slfnResults_48641=PROcapScaled(seqSLFN_48641$V1, seqSLFN_48641$V2, seqSLFN_48641$V5, seqSLFN_48641$V7, extractLoci(loci)[1], extractLoci(loci)[2], bin, windowSize,dist,log(trans))
    fileName=paste(paste(paste(paste("slfn_bin",bin,sep = ""),"_dist",sep = ""),dist,sep = ""),".png",sep = "")
    venn.diagram(
      x = list(which(slfnResults_48646$viterbiPath==1) , which(slfnResults_51644$viterbiPath==1) , 
               which(slfnResults_48641$viterbiPath==1)),
      category.names = c("ID=48646" , "ID=51644" , "ID=48641"),
      filename = fileName,
      output = TRUE ,
      imagetype="png" ,
      height = 640 , 
      width = 640 , 
      resolution = 300,
      compression = "lzw",
      lwd = 2,
      lty = 'blank',
      fill = c('red', 'blue', 'green'),
      cex = .7,
      fontface="bold",
      fontfamily = "sans",
      cat.cex = 0.45,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.fontfamily = "sans")
    
  }
}

loci="17:"
for (bin in binSize){
  for (dist in dists){
    chr17Results_48646=PROcapScaled(seq17_48646$V1, seq17_48646$V2, seq17_48646$V5, seq17_48646$V7, extractLoci(loci)[1], extractLoci(loci)[2], bin, windowSize,dist,log(trans))
    chr17Results_51644=PROcapScaled(seq17_51644$V1, seq17_51644$V2, seq17_51644$V5, seq17_51644$V7, extractLoci(loci)[1], extractLoci(loci)[2], bin, windowSize,dist,log(trans))
    chr17Results_48641=PROcapScaled(seq17_48641$V1, seq17_48641$V2, seq17_48641$V5, seq17_48641$V7, extractLoci(loci)[1], extractLoci(loci)[2], bin, windowSize,dist,log(trans))
    fileName=paste(paste(paste(paste("chr17_bin",bin,sep = ""),"_dist",sep = ""),dist,sep = ""),".png",sep = "")
    venn.diagram(
      x = list(which(chr17Results_48646$viterbiPath==1) , which(chr17Results_51644$viterbiPath==1) , 
               which(chr17Results_48641$viterbiPath==1)),
      category.names = c("ID=48646" , "ID=51644" , "ID=48641"),
      filename = fileName,
      output = TRUE ,
      imagetype="png" ,
      height = 640 , 
      width = 640 , 
      resolution = 300,
      compression = "lzw",
      lwd = 2,
      lty = 'blank',
      fill = c('red', 'blue', 'green'),
      cex = .7,
      fontface="bold",
      fontfamily = "sans",
      cat.cex = 0.45,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      cat.fontfamily = "sans")
    
  }
}

#move away from fixed distance
#1 problem - sliding window too small learn from global distribution genome wide
#use python to summarize data 
#10 bp binning is ok
#do we want to set distance btw reverse and forward
#previously used one peak as anchor and declare peak in window if under certain distance
#if variable peak length need dist
#can classify as single clear peak (each have more than 5 read in 2 kb window), multiple peaks w ambig bds,
#variable peaks (a lot of open chromatin)
#can we categorize using convolution rules
#convolution fn w min value integrate to find max (around 300) learn across genome by sample
#load data efficiently should work
#represent like a marix each position
#w is distance in picture
#concordant btw samples
#convolution nextwork deep learning
#deep bind 3blue1brown