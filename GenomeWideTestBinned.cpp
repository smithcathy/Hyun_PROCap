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
  while(chr(index)==chromosome){
    binIndex=floor(startSite(index)/binSize);
    binnedReads(binIndex,0)+=forwardReads(index);
    binnedReads(binIndex,1)+=revReads(index);
    forLambda+=forwardReads(index);;
    revLambda+=revReads(index);
    index++;
    if (index==size-1){
      break;
    }
  }
}

void OptimDist(IntegerMatrix& rawReads, const int interval, NumericVector& readsPerDistLog, 
               NumericVector& readsPerDistSqrt, const unsigned int maxLen){
  double maxValue=0;
  for(unsigned int i=0; i<maxLen; i++){
    for(int j=-interval; j<interval; j++){
      if (i>abs(j) && i+abs(j)<maxLen){
        if(rawReads(i+j,0)>0 && rawReads(i-j,1)>0){
          readsPerDistLog[j+interval]+=log(rawReads(i+j,0)*rawReads(i-j,1)+1);
          readsPerDistSqrt[j+interval]+=sqrt(rawReads(i+j,0)*rawReads(i-j,1));
        }
      }
    }
  }
}

// [[Rcpp::export]]
List GenomeWideDistBin(IntegerVector& chr, IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads,
                    unsigned int interval, unsigned int maxChrLen, IntegerVector& chrUsed, IntegerVector& maxLenByChr,
                    const int binSize, const int totalLen){
  int maxDist;
  double forLambda, revLambda;
  unsigned int chrCount=chrUsed.size();
  NumericVector readsPerDistLog(interval*2), readsPerDistSqrt(interval*2);
  IntegerMatrix binnedReads(ceil(maxChrLen/binSize),2);
  unsigned int index=0;
  unsigned int size=chr.size();
  for (unsigned int i=0; i<chrCount; i++){
    BinReads(chr,startSite,forwardReads,revReads,binSize,binnedReads,forLambda,revLambda,chrUsed(i),index,size);
    OptimDist(binnedReads, interval, readsPerDistLog, readsPerDistSqrt, ceil(maxLenByChr(i)/binSize));
  }
  double maxValue=0;
  for (unsigned int i=0; i<interval*2; i++){
    if (readsPerDistLog[i]>maxValue){
      maxValue=readsPerDistLog[i];
      maxDist=i-interval;
    }
  }
  forLambda=ceil((double)totalLen/binSize)/forLambda;
  revLambda=ceil((double)totalLen/binSize)/revLambda;
  return List::create(Named("readsPerDistLog")=readsPerDistLog,
                      Named("readsPerDistSqrt")=readsPerDistSqrt,
                      Named("maxDist")=maxDist,
                      Named("forLambda")=forLambda,
                      Named("revLambda")=revLambda);
}