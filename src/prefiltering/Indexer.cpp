#include "Indexer.h"
Indexer::Indexer(const int alphabetSize, const int maxKmerSize){
    this->maxKmerSize = maxKmerSize;
    this->powers = new int[maxKmerSize];
    this->alphabetSize = alphabetSize;
    int pow = 1;
    for( int i=0; i<maxKmerSize; ++i ){
        this->powers[i] = pow;
        pow *= alphabetSize;
    }
    this->maxKmerIndex = 0;
    for( int i=0; i<maxKmerSize; ++i ){
        this->maxKmerIndex += alphabetSize*this->powers[i];
    }

    this->lastKmerIndex = this->maxKmerIndex;

    workspace = new int[100];
}

Indexer::~Indexer(){
    delete[] this->powers;
    delete[] workspace;
}

unsigned int Indexer::int2index( const int *int_seq,const int begin,const int end){
    this->lastKmerIndex = 0;
    for( int i=begin; i<end; i++ ) {
            this->lastKmerIndex += int_seq[i]*this->powers[i-begin];
    }
    return this->lastKmerIndex;
}

unsigned int Indexer::int2index( const int *int_seq){
    int2index(int_seq,0,this->maxKmerSize);
    return this->lastKmerIndex;
}

void Indexer::index2int(int* int_seq, unsigned int idx, int kmerSize){
    for (int i = kmerSize - 1; i >= 0; i--){
        int_seq[i] = idx / powers[i];
        idx = idx - int_seq[i] * powers[i];
    }
}

unsigned int Indexer::getNextKmerIndex (const int* kmer, int kmerSize){
    if (this->lastKmerIndex == this->maxKmerIndex)
        return int2index(kmer, 0, kmerSize);
    else{
        this->lastKmerIndex /= this->alphabetSize;
        this->lastKmerIndex += kmer[kmerSize-1] * this->powers[kmerSize-1];
        return this->lastKmerIndex;
    }
}

void Indexer::reset(){
    this->lastKmerIndex = this->maxKmerIndex;
}

void Indexer::printKmer(int kmerIdx, int kmerSize, char* int2aa){
    index2int(workspace, kmerIdx, kmerSize);
    for (int j = 0; j < kmerSize; j++)
        std::cout << int2aa[workspace[j]];
}

void Indexer::printKmer(const int* kmer, int kmerSize, char* int2aa){
    for (int j = 0; j < kmerSize; j++)
        std::cout << int2aa[kmer[j]];
}

