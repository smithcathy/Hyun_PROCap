#include <Rcpp.h>
#include <fstream>
#include <queue>
using namespace Rcpp;
using namespace std;

void ClearQueue(queue <unsigned int>& full){ //swaps full queue with an empty queue and frees up the memory
  queue <unsigned int> empty;
  swap(full,empty);
}

void OpenFile(string fileName, vector<vector<unsigned int> > &binnedFile, const unsigned int binSize) {
  //cout << "entered OpenFile" << endl;
  string s1;
  unsigned int chr, loci, junk, forward, reverse, index=0;
  vector <unsigned int> row(4);
  ifstream fs(fileName);
  getline(fs, s1);
  while (!getline(fs,s1).eof()){
    stringstream ss(s1);
    ss >> chr >> loci >> junk >> junk >> forward >> junk >> reverse;
    loci=floor(loci/binSize);
    //Rprintf("chr: %d, loci : %d, forward : %d, reverse : %d\n", chr, loci, forward, reverse);
    if (binnedFile.empty()){
      row.push_back(chr);
      row.push_back(loci);
      row.push_back(forward);
      row.push_back(reverse);
      binnedFile.push_back(row);
      row.clear();
    }else{
      //cout << binnedFile.back()[1] << endl;
      //cout << binnedFile.back().front() << endl;
      //cout << index << endl;
      if (binnedFile.back()[1]==loci && binnedFile.back().front()==chr){
        binnedFile[index][2]+=forward;
        binnedFile[index][3]+=reverse;
      }else{
        row.push_back(chr);
        row.push_back(loci);
        row.push_back(forward);
        row.push_back(reverse);
        binnedFile.push_back(row);
        row.clear();
        index++;
      }
    }
  }
  fs.close();
}

void SingleSample(vector<vector<unsigned int> > &binnedFile, vector< vector <double> > &singleSample, 
                  const unsigned int binDist, const unsigned int fileCount){
  //cout << "entered SingleSample" << endl;
  const unsigned int size=binnedFile.size();
  queue <unsigned int> revReadLoci, revReadCount;
  vector <double> row(3);
  unsigned int currentForReadPos, currentForChr, currentForReads, currentRevReads, currentDist, prevChr=0;
  for (unsigned int i=0; i<size; i++){
    currentForChr=binnedFile[i][0];
    currentForReadPos=binnedFile[i][1];
    currentRevReads=binnedFile[i][3];
    if (currentForChr!=prevChr){ //new chromosome - clear queue
      ClearQueue(revReadLoci);
      ClearQueue(revReadCount);
    }
    while (!revReadLoci.empty()){
      currentDist=currentForReadPos-revReadLoci.front();
      if (currentDist>2*binDist){ //reads too far apart - remove rev reads from queue
        revReadLoci.pop();
        revReadCount.pop();
      }else if(currentDist<2*binDist){ //reads too close together - break while loop
        break;
      }else if(currentDist==2*binDist){
        currentForReads=binnedFile[i][2];
        if (currentForReads>0){
          row.push_back(currentForChr);
          row.push_back(currentForReadPos-binDist);
          row.push_back(log(currentForReads*revReadCount.front()+1.0)/fileCount);
          singleSample.push_back(row);
          row.clear();
        }
        revReadLoci.pop();
        revReadCount.pop();
        break;
      }
    }
    prevChr=currentForChr;
    if (currentRevReads>0){
      revReadLoci.push(currentForReadPos);
      revReadCount.push(currentRevReads);
    }
  }
}

void MergeSamples(vector <vector <double> > &singleSample, vector <vector <double> > &mergedSamples){
  //cout << "entered MergeSamples" << endl;
  unsigned int i=0, j=0, sizei=singleSample.size(), sizej=mergedSamples.size();
  double singleReads;
  vector <double> singleRow, mergedRow(5);
  while (i<sizei && j<sizej){
    if (singleRow.empty()){
      singleRow=singleSample[i];
    }
    if (mergedRow.empty()){
      mergedRow=mergedSamples[j];
    }
    if (singleRow.front()==mergedRow.front()){
      if (singleRow[1]<mergedRow[1]){
        singleRow.push_back(singleRow.back());
        singleRow.push_back(1);
        mergedSamples.insert(mergedSamples.begin()+j,singleRow);
        singleRow.clear();
        i++;
      }else if (singleRow[1]>mergedRow[1]){
        mergedRow.clear();
        j++;
      }else{
        singleReads=singleRow.back();
        mergedSamples[j][2]+=singleReads;
        if (mergedSamples[j][3] < singleReads){
          mergedSamples[j][3]=singleReads;
        }
        mergedSamples[j].back()++;
        singleRow.clear();
        mergedRow.clear();
        j++;
        i++;
      }
    }
    else if(singleRow.front()>mergedRow.front()){
      mergedRow.clear();
      j++;
    }
    else if(singleRow.front()<mergedRow.front()){
      singleRow.push_back(singleRow.back());
      singleRow.push_back(1);
      mergedSamples.insert(mergedSamples.begin()+j,singleRow);
      singleRow.clear();
      i++;
    }
  }
  while (i<sizei){
    singleRow=singleSample[i++];
    singleRow.push_back(singleRow.back());
    singleRow.push_back(1);
    mergedSamples.push_back(singleRow);
    singleRow.clear();
  }
}

//void ProcessSamples(vector < vector<double> > &mergedSamples, const unsigned int winDist){
  //unsigned int i, currentChr, currentLoci, size=mergedSamples.size(), prevChr=0;
  //queue <unsigned int> frontWindowLoci, backWindowLoci;
  //double frontSum, frontMax, backSum, backMax;
  //for (i=0; i<winDist; i++){
  //}
  //for(i=0; i<size; i++){
    //currentChr=mergedSamples[i].front();
    //currentLoci=mergedSamples[i][1];
    //if (currentForChr!=prevChr){ //new chromosome - clear queue
      //ClearQueue(revReadLoci);
      //ClearQueue(revReadCount);
    //}
  //}
//}

//[[Rcpp::export]]
vector<vector<double> > MultiSampleHMM(unsigned int binSize, unsigned int windowSize, unsigned int dist, StringVector files){
  unsigned int i, j, size, fileCount=files.size();
  unsigned int binDist=floor(dist/binSize);
  unsigned int winDist=floor(windowSize/binSize);
  vector <double> row;
  vector<vector<unsigned int> > binnedFile;
  vector<vector<double> > mergedSamples, singleSample;
  string file;
  //cout << "intitialized" << endl;
  for (i=0; i<fileCount; i++){
    //cout << i << endl;
    file= as <std::string> (files(i));
    OpenFile(file,binnedFile,binSize);
    SingleSample(binnedFile,singleSample,binDist,fileCount);
    binnedFile.clear();
    if (i>0){
      MergeSamples(singleSample,mergedSamples);
    }else{
      size=singleSample.size();
      for (j=0; j<size; j++){
        row=singleSample.front();
        row.push_back(row.back());
        row.push_back(1);
        mergedSamples.push_back(row);
        row.clear();
        singleSample.erase(singleSample.begin());
      }
    }
    singleSample.clear();
  }
  //ProcessSamples(mergedSamples,winDist);
  return mergedSamples;
}