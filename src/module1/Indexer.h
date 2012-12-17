#ifndef INDEXER_H
#define INDEXER_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de, Martin Steinegger Martin.Steinegger@campus.lmu.de
// 
// Manages the conversion of the int coded k-mer into a int index and vice versa.
//


#include <string>
#include <iostream>

class Indexer{
	
	public:
        Indexer(const size_t alphabetSize, const size_t maxKmerSize);
        ~Indexer();

        // get the index of the k-mer, beginning at "begin" in the int_seq and ending at "end"
        size_t int2index( const int *int_seq,const int begin,const int end);
        // get the index of the k-mer of length maxKmerSize, beginning at position 0
        size_t int2index( const int *int_seq);
        
        // get the int sequence for the k-mer with the index idx of kmerSize
        void index2int(int* int_seq, int idx, int kmerSize);
        
        size_t getNextKmerIndex(int* kmer, int kmerSize);
        
        size_t * powers;

    private:

        size_t alphabetSize;
        
        size_t maxKmerSize;

        int lastKmerIndex;
};
#endif
