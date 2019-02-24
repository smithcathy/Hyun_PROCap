#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

void BinReads(IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads, unsigned int lociStart, 
              unsigned int lociLen, unsigned int binSize, IntegerMatrix& binnedReads){
  unsigned int size=(int)startSite.size();
  unsigned int binNumber=binnedReads.nrow();
  unsigned int m=0;
  //forLambda, revLambda=0;
  fill(binnedReads.begin(),binnedReads.end(),0);
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
    //forLambda+=binnedReads(i,0);
    //revLambda+=binnedReads(i,1);
  }
  //forLambda=1/(forLambda/binNumber);
  //revLambda=1/(revLambda/binNumber);
}


double GetLambda(int* winStart, unsigned int start, unsigned int end){
  double sum=0.0;
  for (unsigned int i=start; i<end; i++){
    sum+=*winStart;
    winStart++;
  }
  if (sum==0){
    return 0; //avoids division by 0
  } else{
    return 1/(sum/(end-start)); 
  }
}

void ObsProb(IntegerMatrix& binnedReads, unsigned int binSize, unsigned int windowSize, NumericMatrix& obsProb, unsigned int binDist){
  unsigned int binsPerWindow=windowSize/binSize;
  unsigned int size=binnedReads.nrow()-binDist; //truncates edges that cannot contain state 1 a priori
  if (binDist>size){
    stop("Sequence must be longer than double the specified distance");
  }
  if (size+binDist<binsPerWindow){
    stop("Sequence must be longer than the specified window size");
  }
  double p=log(.001); //lower bd for probabilities
  double q=log(.99); //upper bd for probabilities
  int* startFor; 
  int* startRev;
  double lambdaFor, lambdaRev, obsForS0, obsForS1, obsRevS0, obsRevS1;
  for (unsigned int i=0; i<size; i++){ // when rev reads and for reads are in the middle and can take the lambda normally
    if (i>=binsPerWindow/2.0 && i+binsPerWindow/2.0<=size){
      startFor=&binnedReads(i+binDist-binsPerWindow/2,0);
      startRev=&binnedReads(i-binsPerWindow/2,1);
      lambdaFor=GetLambda(startFor,i+binDist-binsPerWindow/2,i+binDist+binsPerWindow/2);
      lambdaRev=GetLambda(startRev,i-binsPerWindow/2,i+binsPerWindow/2);
    }else if (i<binsPerWindow/2.0){ //not enough bins behind i for full window
      startRev=&binnedReads(0,1);
      lambdaRev=GetLambda(startRev,0,i+binsPerWindow/2);
      if (i+binDist>=binsPerWindow/2.0){ //enough bins behind for full forward window
        startFor=&binnedReads(i+binDist-binsPerWindow/2,0);
        lambdaFor=GetLambda(startFor,i+binDist-binsPerWindow/2,i+binDist+binsPerWindow/2);
      }else{ //not enough bins behind for forward window
        startFor=&binnedReads(0,0);
        lambdaFor=GetLambda(startFor,0,i+binDist+binsPerWindow/2);
      }
    }else{ //not enough bins ahead for full window
      startFor=&binnedReads(i+binDist-binsPerWindow/2,0);
      lambdaFor=GetLambda(startFor,i+binDist-binsPerWindow/2,size+binDist);
      if (i+windowSize/2.0<=size+binDist){ //enough reads in front of rev for full window
        startRev=&binnedReads(i-binsPerWindow/2,1);
        lambdaRev=GetLambda(startRev,i-binsPerWindow/2,i+binsPerWindow/2);
      }else{ //not enough reads in front of rev for full window
        startRev=&binnedReads(i-binsPerWindow/2,1);
        lambdaRev=GetLambda(startRev,i-binsPerWindow/2,size+binDist);
      }
    }
    obsForS1=min(max(log(1-exp(-lambdaFor*binnedReads(i+binDist,0))),p),q); //if prob is 0 set to p instead
    obsRevS1=min(max(log(1-exp(-lambdaRev*binnedReads(i,1))),p),q);
    obsForS0=min(max(-.25*binnedReads(i+binDist,0),p),q); // inverse cdf (survival): log(1-(1-exp(-lambda*binnedReads)))==-lambda*binnedReads
    obsRevS0=min(max(-.25*binnedReads(i,1),p),q); //keep q>=prob>=p
    obsProb(i,0)=obsForS0+obsRevS0;
    obsProb(i,1)=obsForS1+obsRevS1;
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
    delta(i,0)=pi(i)+obsProb(0,i);
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
      delta(i,t)+=obsProb(t,i);
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
List PROcapScaled(IntegerVector chr, IntegerVector startSite, IntegerVector forwardReads, IntegerVector revReads, unsigned int lociStart, 
            unsigned int lociLen, unsigned int binSize, unsigned int windowSize, unsigned int dist, NumericMatrix trans){
  // if binSize not divisible by lociLen, add one bin to cover remaining bases
  if (trans.nrow()!=trans.ncol()){
    stop("Transition matrix must be square");
  }
  unsigned int states=trans.nrow();
  unsigned int binNumber=ceil(lociLen/binSize);
  unsigned int binDist=ceil(dist/binSize);
  if(binSize>dist){
    stop("Distance between TFBS must be greater than binning length");
  }
  IntegerMatrix binnedReads(binNumber,states);
  NumericMatrix obsProb(binNumber-binDist,states), condProb(binNumber-binDist,states);
  IntegerVector viterbiPath(binNumber-binDist);
  BinReads(startSite,forwardReads,revReads,lociStart,lociLen,binSize,binnedReads);
  //print(binnedReads);
  ObsProb(binnedReads, binSize, windowSize, obsProb, binDist);
  //print(obsProb);
  NumericVector pi=GetPi(trans);
  ViterbiPath(obsProb, viterbiPath, trans, pi);
  //print(viterbiPath);
  for (unsigned int i=0; i<states; ++i){
    obsProb(_,i)=exp(obsProb(_,i));
    trans(_,i)=exp(trans(_,i));
  }
  pi=exp(pi);
  ForwardBackward(obsProb, condProb, trans, pi);
  return List::create(Named("binnedReads")=binnedReads,
                       Named("obsProb")=obsProb,
                       Named("viterbiPath")=viterbiPath,
                       Named("condProb")=condProb);
}