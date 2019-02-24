#include <Rcpp.h>
#include <iostream>
#include <vector>
using namespace Rcpp;
using namespace std;

void BinReads(IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads, const int lociStart, 
              const int lociLen, const int binSize, IntegerMatrix& binnedReads, double& forLambda,
              double& revLambda){
  unsigned int size=(int)startSite.size();
  unsigned int binNumber=binnedReads.nrow();
  unsigned int m=0;
  for (unsigned int i=0; i<binNumber; i++){
    unsigned int start=lociStart+i*binSize;
    while(m<size){
      if (startSite(m)>=start+binSize){  //deals with case where no reads in bin
        break; //and breaks out of while loop without incrementing
      } else if (startSite(m)>=start & startSite(m)<start+binSize){
        binnedReads(i,0)+=forwardReads(m); //first column is forward reads
        binnedReads(i,1)+=revReads(m); //second column is backward reads
        m++;
      }
    }
    forLambda+=binnedReads(i,0);
    revLambda+=binnedReads(i,1);
  }
  forLambda=binNumber/forLambda;
  revLambda=binNumber/revLambda;
}

void CalcProb(IntegerMatrix& binnedReads,NumericMatrix& prob,const double& forLambda,
              const double& revLambda){
  unsigned int bins=binnedReads.nrow();
  double p=.001; //lower bd for probabilities
  double q=.99; //upper bd for probabilities
  for (unsigned int i=0; i<bins; i++){
    prob(i,0)=min(max(1-exp(-forLambda*binnedReads(i,0)),p),q);
    prob(i,1)=min(max(1-exp(-revLambda*binnedReads(i,1)),p),q);
    prob(i,2)=min(max(exp(-.25*binnedReads(i,0)),p),q); // inverse cdf (survival): log(1-(1-exp(-lambda*binnedReads)))==-lambda*binnedReads
    prob(i,3)=min(max(exp(-.25*binnedReads(i,1)),p),q);
  }
}

void ConvolDist(IntegerMatrix& binnedReads, NumericMatrix& prob, const int binInterval, NumericVector& readsPerDist,
                unsigned int& maxDist){
  unsigned int bins=binnedReads.nrow();
  double maxValue=-1;
  for (int i=0; i<bins; i++){
    for (int j=-binInterval; j<binInterval; j++){
      if (i>=-j && i+j<bins){
        readsPerDist(j+binInterval)+=prob(i,0)*prob(i+j,1);
        if (i==bins-1 && readsPerDist(j+binInterval)>maxValue){
          maxDist=j+binInterval;
          maxValue=readsPerDist(j+binInterval);
        }
      }
    }
  }
}

void ObsProb(IntegerMatrix& binnedReads, NumericMatrix& prob, unsigned int binSize, NumericMatrix& obsProb, const int maxDist){
  unsigned int size=obsProb.nrow();//truncates edges that cannot contain state 1 a priori
  for (unsigned int i=0; i<size; i++){
    obsProb(i,0)=prob(i+maxDist,2)*prob(i,3);
    obsProb(i,1)=prob(i+maxDist,0)*prob(i,1);
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
List PROcapBestDist(IntegerVector chr, IntegerVector startSite, IntegerVector forwardReads, IntegerVector revReads, unsigned int lociStart, 
                unsigned int lociLen, unsigned int binSize, unsigned int interval, NumericMatrix trans){
  if (trans.nrow()!=trans.ncol()){
    stop("Transition matrix must be square");
  }
  const int states=trans.nrow();
  // if binSize not divisible by lociLen, add one bin to cover remaining bases
  const int binNumber=ceil(lociLen/binSize);
  const int binInterval=ceil(interval/binSize);
  unsigned int maxDist;
  double forLambda, revLambda;
  NumericVector readsPerDist(binInterval*2);
  NumericMatrix prob(binNumber,2*states);
  IntegerMatrix binnedReads(binNumber,states);
  BinReads(startSite,forwardReads,revReads,lociStart,lociLen,binSize,binnedReads,forLambda,revLambda);
  CalcProb(binnedReads,prob,forLambda,revLambda);
  ConvolDist(binnedReads,prob,binInterval,readsPerDist, maxDist);
  maxDist=binInterval-maxDist;
  NumericMatrix obsProb(binNumber-maxDist,states), condProb(binNumber-maxDist,states);
  IntegerVector viterbiPath(binNumber-maxDist);
  ObsProb(binnedReads, prob, binSize, obsProb, maxDist);
  NumericVector pi=GetPi(trans);
  ViterbiPath(obsProb, viterbiPath, trans, pi);
  pi=exp(pi);
  ForwardBackward(obsProb, condProb, trans, pi);
  return List::create(Named("binnedReads")=binnedReads,
                      Named("readsPerDist")=readsPerDist,
                      Named("maxDist")=maxDist,
                      Named("forLambda")=forLambda,
                      Named("revLambda")=revLambda,
                      Named("obsProb")=obsProb,
                      Named("viterbiPath")=viterbiPath,
                      Named("condProb")=condProb);
}

