#include <Rcpp.h>
#include <iostream>
#include <vector>
using namespace Rcpp;
using namespace std;

void BinReads(IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads, const int lociStart, 
              const int lociLen, const int binSize, IntegerMatrix& binnedReads, IntegerMatrix& rawReads, 
              double& forLambda, double& revLambda){
  unsigned int size=(int)startSite.size();
  unsigned int binNumber=binnedReads.nrow();
  unsigned int rawIndex, binIndex;
  for (unsigned int i=0; i<size; i++){
    rawIndex=startSite(i)-lociStart;
    binIndex=floor((startSite(i)-lociStart)/binSize);
    rawReads(rawIndex,0)=forwardReads(i);
    rawReads(rawIndex,1)=revReads(i);
    binnedReads(binIndex,0)+=forwardReads(i);
    binnedReads(binIndex,1)+=revReads(i);
    forLambda+=forwardReads(i);;
    revLambda+=revReads(i);
  }
  forLambda=binNumber/forLambda;
  revLambda=binNumber/revLambda;
}

void OptimDist(IntegerMatrix& rawReads, const int interval, NumericVector& readsPerDistLog, 
               NumericVector& readsPerDistSqrt, int& maxDist){
  unsigned int size=rawReads.nrow();
  double maxValue=0;
  for(unsigned int i=0; i<size; i++){
    for(int j=-interval; j<interval; j++){
      if (i==(size-1) && readsPerDistLog[j+interval]>maxValue){
        maxValue=readsPerDistLog[j+interval];
        maxDist=j+1;
      }
      if (i>abs(j) && i+abs(j)<size){
        if(rawReads(i+j,0)>0 && rawReads(i-j,1)>0){
          readsPerDistLog[j+interval]+=log(rawReads(i+j,0)*rawReads(i-j,1)+1);
          readsPerDistSqrt[j+interval]+=sqrt(rawReads(i+j,0)*rawReads(i-j,1));
        }
      }
    }
  }
}

void ObsProb(IntegerMatrix& binnedReads, NumericMatrix& obsProb, const int distBin, const double& forLambda, 
             const double& revLambda){
  double p=.0001; //lower bd for probabilities
  double q=.999; //upper bd for probabilities
  unsigned int size=obsProb.nrow();//truncates edges that cannot contain state 1 a priori
  if (distBin>0){
    for (unsigned int i=0; i<size; i++){
      obsProb(i,0)=min(max(.6025*exp(binnedReads(i+2*distBin,0)+binnedReads(i,1)),p),q);
      obsProb(i,1)=min(max(1-exp(-forLambda*binnedReads(i+2*distBin,0))*(1-exp(-revLambda*binnedReads(i,1))),p),q);
    }
  }else if (distBin<=0){
    for (unsigned int i=0; i<size; i++){
      obsProb(i,0)=min(max(.6025*exp(binnedReads(i,0)+binnedReads(i+2*distBin,1)),p),q);
      obsProb(i,1)=min(max(1-exp(-forLambda*binnedReads(i,0))*(1-exp(-revLambda*binnedReads(i+2*distBin,1))),p),q);
    }
  }
  
}

void ObsProb(IntegerMatrix& binnedReads, NumericMatrix& prob, unsigned int binSize, NumericMatrix& obsProb, const int distBin){
  unsigned int size=obsProb.nrow();//truncates edges that cannot contain state 1 a priori
  for (unsigned int i=0; i<size; i++){
    obsProb(i,0)=prob(i+distBin,2)*prob(i,3);
    obsProb(i,1)=prob(i+distBin,0)*prob(i,1);
  }
}


NumericVector GetPi(NumericMatrix& trans){
  unsigned int states=trans.nrow();
  NumericVector startPi(states), newPi(states);
  if (trans.nrow()!=trans.ncol()){
    stop("The transition matrix should be square");
  }
  for (unsigned int i=0; i<states; i++){
    startPi(i)=1.0/states;
  }
  for (unsigned int i=0; i<states; i++){
    for (unsigned int j=0; j<states; j++){
      newPi(i)+=exp(trans(i,j))*startPi(i);
    }
  }
  return log(newPi);
}

void ViterbiPath(NumericMatrix& obsProb, IntegerVector& viterbiPath, NumericMatrix& trans,
                 NumericVector& pi) {
  unsigned int n = obsProb.nrow();
  unsigned int states=obsProb.ncol();
  NumericMatrix delta(states,n);
  IntegerMatrix phi(states,n);
  double maxj=R_NegInf;
  for (unsigned int i=0; i<states; i++){
    delta(i,0)=pi(i)+log(obsProb(0,i));
  }
  for (unsigned int t=1; t<n; ++t){
    for (unsigned int i=0; i<states; ++i){
      maxj=R_NegInf;
      for (unsigned int j=0; j<states; ++j){
        double v=delta(j,t-1)+trans(j,i);
        if (v>maxj){
          maxj=v;
          delta(i,t)=v;
          phi(i,t)=j;
        }
      }
      delta(i,t)+=log(obsProb(t,i));
    }
  }
  double max=R_NegInf;
  for (unsigned int i=0; i<states; i++){
    if (delta(i,n-1)>max){
      max=delta(i,n-1);
      viterbiPath(n-1)=i;
    }
  }
  for (unsigned int j=n-1; j>0; --j){
    viterbiPath[j-1]=phi(viterbiPath[j],j);
  }
}

void ForwardBackward(NumericMatrix& obsProb, NumericMatrix& condProb, NumericMatrix& trans, NumericVector& pi){
  unsigned int n=obsProb.nrow();
  unsigned int states=pi.size();
  NumericMatrix alpha(states,n), beta(states,n);
  NumericVector scale(n);
  for (unsigned int i=0; i<states; ++i){
    scale(0)+=pi(i)*obsProb(0,i);
  }
  for (unsigned int i=0; i<states; ++i){
    alpha(i,0)=pi(i)*obsProb(0,i)/scale(0);
    beta(i,n-1)=1;
  }
  for (unsigned int t=1; t<n; ++t){
    scale(t)=0;
    for (unsigned int i=0; i<states; ++i){
      for (unsigned int j=0; j<states; ++j){
        scale(t)+=alpha(j,t-1)*trans(j,i);
      }
      scale(t)*=obsProb(t,i);
    }
    for (unsigned int i=0; i<states; ++i){
      alpha(i,t)=0;
      for (unsigned int j=0; j<states; ++j){
        alpha(i,t)+=alpha(j,t-1)*trans(j,i)/scale(t);
      }
      alpha(i,t)*=obsProb(t,i);
    }
  }
  for (int t=n-2; t>=0; --t){
    for (unsigned int i=0; i<states; ++i){
      beta(i,t)=0;
      for (unsigned int j=0; j<states; ++j){
        beta(i,t)+=beta(j,t+1)*trans(i,j)*obsProb(t+1,j)/scale(t+1);
      }
    }
  }
  for (unsigned int t=0; t<n; ++t){
    double sum=0;
    for (unsigned int i=0; i<states; ++i){
      sum+=alpha(i,t)*beta(i,t); 
    }
    for (unsigned int i=0; i<states; ++i){
      condProb(t,i)=(alpha(i,t)*beta(i,t))/sum; 
    }
  }
}

// [[Rcpp::export]]
List PROcapOptimDist(IntegerVector chr, IntegerVector startSite, IntegerVector forwardReads, IntegerVector revReads, unsigned int lociStart, 
                    unsigned int lociLen, unsigned int binSize, unsigned int interval, NumericMatrix trans){
  if (trans.nrow()!=trans.ncol()){
    stop("Transition matrix must be square");
  }
  const int states=trans.nrow();
  // if binSize not divisible by lociLen, add one bin to cover remaining bases
  const int binNumber=ceil(lociLen/binSize);
  //const int binInterval=ceil(interval/binSize);
  int maxDist;
  double forLambda, revLambda;
  NumericVector readsPerDistLog(interval*2), readsPerDistSqrt(interval*2);
  IntegerMatrix binnedReads(binNumber,states), rawReads(lociLen,states);
  BinReads(startSite,forwardReads,revReads,lociStart,lociLen,binSize,binnedReads,rawReads,forLambda,revLambda);
  OptimDist(rawReads, interval, readsPerDistLog, readsPerDistSqrt, maxDist);
  unsigned int distBin=ceil(maxDist/binSize);
  NumericMatrix obsProb(binNumber-2*distBin,states), condProb(binNumber-2*distBin,states);
  IntegerVector viterbiPath(binNumber-2*distBin);
  ObsProb(binnedReads, obsProb, distBin, forLambda,revLambda);
  NumericVector pi=GetPi(trans);
  ViterbiPath(obsProb, viterbiPath, trans, pi);
  pi=exp(pi);
  ForwardBackward(obsProb, condProb, trans, pi);
  return List::create(Named("binnedReads")=binnedReads,
                      Named("readsPerDistLog")=readsPerDistLog,
                      Named("readsPerDistSqrt")=readsPerDistSqrt,
                      Named("maxDist")=maxDist,
                      Named("forLambda")=forLambda,
                      Named("revLambda")=revLambda,
                      Named("obsProb")=obsProb,
                      Named("viterbiPath")=viterbiPath,
                      Named("condProb")=condProb);
}