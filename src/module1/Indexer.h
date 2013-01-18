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
        Indexer(const int alphabetSize, const int maxKmerSize);
        ~Indexer();

        // get the index of the k-mer, beginning at "begin" in the int_seq and ending at "end"
        int int2index( const int *int_seq,const int begin,const int end);
        // get the index of the k-mer of length maxKmerSize, beginning at position 0
        int int2index( const int *int_seq);
        
        // get the int sequence for the k-mer with the index idx of kmerSize
        void index2int(int* int_seq, int idx, int kmerSize);
        
        int getNextKmerIndex(int* kmer, int kmerSize);
        
        int * powers;

    private:

        int alphabetSize;
        
        int maxKmerSize;

        int lastKmerIndex;
};
#endif
