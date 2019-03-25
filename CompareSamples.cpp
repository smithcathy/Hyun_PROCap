#include <Rcpp.h>
#include <queue>
using namespace Rcpp;
using namespace std;

void ClearQueue(queue <unsigned int>& full){ //swaps full queue with an empty queue and frees up the memory
  queue <unsigned int> empty;
  swap(full,empty);
}

void ReadsatDistance(IntegerVector& chr, IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads, 
                     const int dist, NumericMatrix& reads, double& forLambda, double& revLambda){
  const unsigned int size=chr.size();
  queue <unsigned int> revReadLoci, revReadCount;
  unsigned int currentForReadPos, currentDist, readsIndex=0, prevChr=0;
  for (unsigned int i=0; i<size; i++){
    currentForReadPos=startSite(i);
    if (chr(i)!=prevChr){ //new chromosome - clear queue
      ClearQueue(revReadLoci);
      ClearQueue(revReadCount);
    }
    while (!revReadLoci.empty()){
      currentDist=currentForReadPos-revReadLoci.front();
      if (currentDist>2*dist){ //reads too far apart - remove rev reads from queue
        revReadLoci.pop();
        revReadCount.pop();
      }else if(currentDist<2*dist){ //reads too close together - break while loop
        break;
      }else if(currentDist==2*dist){
        if (forwardReads(i)>0 && revReadCount.front()>0){
          reads(readsIndex,0)=chr(i);
          reads(readsIndex,1)=currentForReadPos-dist;
          reads(readsIndex,2)=log(forwardReads(i)*revReadCount.front()+1);
          readsIndex++;
        }
        revReadLoci.pop();
        revReadCount.pop();
        break;
      }
    }
    prevChr=chr(i);
    revReadLoci.push(startSite(i));
    revReadCount.push(revReads(i));
    forLambda+=forwardReads(i);
    revLambda+=revReads(i);
  }
}

// [[Rcpp::export]]
List GenomeWideCompare(IntegerVector& chr, IntegerVector& startSite, IntegerVector& forwardReads, IntegerVector& revReads,
                    unsigned int dist, const unsigned int totalLen){
  double forLambda=0, revLambda=0;
  unsigned int size=chr.size();
  NumericMatrix reads(size,3);
  ReadsatDistance(chr, startSite, forwardReads, revReads, dist, reads, forLambda, revLambda);
  forLambda=(double)totalLen/forLambda;
  revLambda=(double)totalLen/revLambda;
  return List::create(Named("reads")=reads,
                      Named("forLambda")=forLambda,
                      Named("revLambda")=revLambda);
}