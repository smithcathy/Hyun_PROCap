#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <map>
using namespace Rcpp;
using namespace std;

void RawReads(IntegerVector& chr, IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads, 
              IntegerMatrix& rawReads, double& forLambda, double& revLambda, const unsigned int chromosome, 
              unsigned int& index, const unsigned int size){
  fill(rawReads.begin(),rawReads.end(),0);
  while (chr(index)==chromosome){
    rawReads(startSite(index),0)=forwardReads(index);
    rawReads(startSite(index),1)=revReads(index);
    forLambda+=forwardReads(index);
    revLambda+=revReads(index);
    index++;
    if (index==(size-1)){
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
List GenomeWideDist(IntegerVector& chr, IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads,
                    unsigned int interval, unsigned int maxChrLen, IntegerVector& chrUsed, IntegerVector& maxLenByChr,
                    const unsigned int totalLen){
  int maxDist;
  double forLambda, revLambda;
  unsigned int chrCount=chrUsed.size();
  NumericVector readsPerDistLog(interval*2), readsPerDistSqrt(interval*2);
  IntegerMatrix rawReads(maxChrLen,2);
  unsigned int index=0;
  unsigned int size=chr.size();
  for (unsigned int i=0; i<chrCount; i++){
    RawReads(chr,startSite,forwardReads,revReads,rawReads,forLambda,revLambda,chrUsed(i),index,size);
    OptimDist(rawReads, interval, readsPerDistLog, readsPerDistSqrt, maxLenByChr(i));
    cout << i << endl;
  }
  double maxValue=0;
  for (unsigned int i=0; i<interval*2; i++){
    if (readsPerDistLog[i]>maxValue){
      maxValue=readsPerDistLog[i];
      maxDist=i-interval;
    }
  }
  forLambda=totalLen/forLambda;
  revLambda=totalLen/revLambda;
  return List::create(Named("readsPerDistLog")=readsPerDistLog,
                      Named("readsPerDistSqrt")=readsPerDistSqrt,
                      Named("maxDist")=maxDist,
                      Named("forLambda")=forLambda,
                      Named("revLambda")=revLambda);
}