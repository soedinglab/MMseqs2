#include "Indexer.h"
#include "Debug.h"

Indexer::Indexer(const size_t alphabetSize, const int maxKmerSize){
    this->maxKmerSize = maxKmerSize;
    this->powers = new size_t[maxKmerSize];
    this->alphabetSize = alphabetSize;
    size_t pow = 1;
    for( int i=0; i<maxKmerSize; ++i ){
        this->powers[i] = pow;
        pow *= alphabetSize;
    }
    this->maxKmerIndex = 0;
    for( int i=0; i<maxKmerSize; ++i ){
        this->maxKmerIndex += alphabetSize*this->powers[i];
    }

    this->lastKmerIndex = this->maxKmerIndex;

    workspace = new size_t[100];
}

Indexer::~Indexer(){
    delete[] this->powers;
    delete[] workspace;
}

size_t Indexer::int2index( const int *int_seq,const int begin,const int end){
    this->lastKmerIndex = 0;
    for(int i = begin; i < end; i++) {
        this->lastKmerIndex += int_seq[i]*this->powers[i-begin];
    }
    return this->lastKmerIndex;
}

size_t Indexer::int2index( const int *int_seq){
    int2index(int_seq, 0, this->maxKmerSize);
    return this->lastKmerIndex;
}


void Indexer::reset(){
    this->lastKmerIndex = this->maxKmerIndex;
}

void Indexer::printKmer(size_t kmerIdx, int kmerSize, char* int2aa){
    index2int(workspace, kmerIdx, kmerSize);
    for (int j = 0; j < kmerSize; j++)
        Debug(Debug::INFO) << int2aa[workspace[j]];
}

void Indexer::printKmer(const int* kmer, int kmerSize, char* int2aa){
    for (int j = 0; j < kmerSize; j++)
        Debug(Debug::INFO) << int2aa[kmer[j]];
}