#include "Indexer.h"
Indexer::Indexer(const size_t alphabetSize, const size_t maxKmerSize){
    this->maxKmerSize = maxKmerSize;
    this->powers = new size_t[maxKmerSize];
    this->alphabetSize = alphabetSize;
    size_t pow = 1;
    for( size_t i=0; i<maxKmerSize; ++i ){
        this->powers[i] = pow;
        pow *= alphabetSize;
    }
    this->lastKmerIndex = -1;
}

Indexer::~Indexer(){
    delete this->powers;
}

size_t Indexer::int2index( const int *int_seq,const int begin,const int end){
    this->lastKmerIndex = 0;
    for( size_t i=begin; i<end; i++ ) {
            this->lastKmerIndex += int_seq[i]*this->powers[i-begin];
    }
    return this->lastKmerIndex;
}

size_t Indexer::int2index( const int *int_seq){
    int2index(int_seq,0,this->maxKmerSize);
    return this->lastKmerIndex;
}

void Indexer::index2int(int* int_seq, int idx, int kmerSize){
    for (int i = kmerSize - 1; i >= 0; i--){
        int_seq[i] = idx / powers[i];
        idx = idx - int_seq[i] * powers[i];
    }
}

size_t Indexer::getNextKmerIndex (int* kmer, int kmerSize){
    if (this->lastKmerIndex == -1)
        return int2index(kmer, 0, kmerSize);
    else{
        this->lastKmerIndex /= this->alphabetSize;
        this->lastKmerIndex += kmer[kmerSize-1] * this->powers[kmerSize-1];
        return this->lastKmerIndex;
    }
}
