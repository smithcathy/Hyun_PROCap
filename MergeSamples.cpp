#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

void MergeandCount(NumericMatrix& sample1, NumericMatrix& sample2, NumericMatrix& combined, IntegerVector& counts){
  unsigned int size1=sample1.nrow();
  unsigned int size2=sample2.nrow();
  unsigned int chri, chrj, locii, locij, i=0, j=0, index=0;
  while (i<size1 || j<size2){
    if(i<size1 && j<size2){
      chri=sample1(i,0);
      chrj=sample2(j,0);
      locii=sample1(i,1);
      locij=sample2(j,1);
      if (chri==chrj){
        if (locii<locij){
          combined(index,0)=chri;
          combined(index,1)=locii;
          combined(index,2)=sample1(i,2);
          combined(index,3)=0;
          counts(1)++;
          i++;
        }else if (locii>locij){
          combined(index,0)=chrj;
          combined(index,1)=locij;
          combined(index,2)=0;
          combined(index,3)=sample2(j,2);
          counts(2)++;
          j++;
        }else{
          combined(index,0)=chri;
          combined(index,1)=locii;
          combined(index,2)=sample1(i,2);
          combined(index,3)=sample2(j,2);
          counts(3)++;
          i++;
          j++;
        }
      }else if (chri<chrj){
        combined(index,0)=chri;
        combined(index,1)=locii;
        combined(index,2)=sample1(i,2);
        combined(index,3)=0;
        counts(1)++;
        i++;
      }else{
        combined(index,0)=chrj;
        combined(index,1)=locij;
        combined(index,2)=0;
        combined(index,3)=sample2(j,2);
        counts(2)++;
        j++;
      }
    }else if (i<size1){
      combined(index,0)=sample1(i,0);
      combined(index,1)=sample1(i,1);
      combined(index,2)=sample1(i,2);
      combined(index,3)=0;
      counts(1)++;
      i++;
    }else{
      combined(index,0)=sample2(j,0);
      combined(index,1)=sample2(j,1);
      combined(index,2)=0;
      combined(index,3)=sample2(j,2);
      counts(2)++;
      j++;
    }
    index++;
  }
  counts(0)=index;
}

// [[Rcpp::export]]
List MergeSamples(NumericMatrix& sample1, NumericMatrix& sample2){
  NumericMatrix combined(sample1.nrow()+sample2.nrow(),4);
  IntegerVector counts(4);
  MergeandCount(sample1, sample2, combined, counts);
  return List::create(Named("mergedData")=combined,
                      Named("counts")=counts);
}