#ifndef INDEXER_H
#define INDEXER_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de, Martin Steinegger martin.steinegger@mpibpc.mpg.de
//
// Manages the conversion of the int coded k-mer into a int index and vice versa.
//


#include <string>
#include <iostream>

class Indexer{
    
public:
    Indexer(const size_t alphabetSize, const int maxKmerSize);
    ~Indexer();
    
    // get the index of the k-mer, beginning at "begin" in the int_seq and ending at "end"
    size_t int2index( const int *int_seq,const int begin,const int end);
    // get the index of the k-mer of length maxKmerSize, beginning at position 0
    size_t int2index( const int *int_seq);
    
    // get the int sequence for the k-mer with the index idx of kmerSize
    inline void index2int(size_t * int_seq, size_t idx, int kmerSize){
        for (int i = kmerSize - 1; i >= 0; i--){
            int_seq[i] = idx / powers[i];
            idx = idx - int_seq[i] * powers[i];
        }
    }
    
    // k-mer iterator, remembers the last k-mer
    size_t getNextKmerIndex (const int* kmer, int kmerSize){
        if (this->lastKmerIndex == this->maxKmerIndex)
            return int2index(kmer, 0, kmerSize);
        else{
            this->lastKmerIndex /= this->alphabetSize;
            this->lastKmerIndex += kmer[kmerSize-1] * this->powers[kmerSize-1];
            return this->lastKmerIndex;
        }
    }
    
    // reset the last k-mer
    void reset();
    
    // print k amino acids of the k-mer with index kmerIdx
    // int k-mer is written into workspace
    void printKmer(size_t kmerIdx, int kmerSize, char* int2aa);
    
    // print k amino acids of int k-mer kmer
    void printKmer(const int* kmer, int kmerSize, char* int2aa);
    
    size_t * powers;
    
private:
    
    size_t alphabetSize;
    
    size_t maxKmerSize;
    
    size_t lastKmerIndex;
    
    size_t maxKmerIndex;
    
    size_t * workspace;
};
#endif
