#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <map>
using namespace Rcpp;
using namespace std;

void BinReads(IntegerVector& chr, IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads, 
              const int binSize, IntegerMatrix& binnedReads, double& forLambda, double& revLambda, 
              const unsigned int chromosome,unsigned int& index, const unsigned int size){
  unsigned int binNumber=binnedReads.nrow();
  fill(binnedReads.begin(),binnedReads.end(),0);
  unsigned int binIndex;
  while(index<size && chr(index)==chromosome){
    binIndex=floor(startSite(index)/binSize);
    binnedReads(binIndex,0)+=forwardReads(index);
    binnedReads(binIndex,1)+=revReads(index);
    forLambda+=forwardReads(index);;
    revLambda+=revReads(index);
    index++;
  }
}

void ReadsAtDist(IntegerMatrix& binnedReads, const int dist, NumericVector& reads, const unsigned int maxLen, 
                 const unsigned int col){
  for (unsigned int i=0; i<maxLen; i++){
    if(binnedReads(i+2*dist,0)>0 && binnedReads(i,1)>0){
      reads(i, col)=log(binnedReads(i+2*dist,0)*binnedReads(i,1)+1);
    }
  }
}

// [[Rcpp::export]]
List GenomeWideDistBin(IntegerVector& chr, IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads,
                    unsigned int dist, unsigned int maxChrLen, IntegerVector& chrUsed, IntegerVector& maxLenByChr,
                    const int binSize, const int totalLen){
  double forLambda=0, revLambda=0;
  unsigned int chrCount=chrUsed.size();
  NumericMatrix reads(ceil(maxChrLen/binSize)-2*dist,chrCount);
  IntegerMatrix binnedReads(ceil(maxChrLen/binSize),2);
  unsigned int index=0;
  unsigned int size=chr.size();
  for (unsigned int i=0; i<chrCount; i++){
    BinReads(chr,startSite,forwardReads,revReads,binSize,binnedReads,forLambda,revLambda,chrUsed(i),index,size);
    ReadsAtDist(binnedReads, dist, reads, ceil(maxLenByChr(i)/binSize)-2*dist, i);
  }
  forLambda=ceil((double)totalLen/binSize)/forLambda;
  revLambda=ceil((double)totalLen/binSize)/revLambda;
  return List::create(Named("reads")=reads,
                      Named("forLambda")=forLambda,
                      Named("revLambda")=revLambda);
}