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
        unsigned int int2index( const int *int_seq,const int begin,const int end);
        // get the index of the k-mer of length maxKmerSize, beginning at position 0
        unsigned int int2index( const int *int_seq);
        
        // get the int sequence for the k-mer with the index idx of kmerSize
        void index2int(int* int_seq, unsigned int idx, int kmerSize);
       
        // k-mer iterator, remembers the last k-mer
        unsigned int getNextKmerIndex(const int* kmer, int kmerSize);

        // reset the last k-mer
        void reset();

        // print k amino acids of the k-mer with index kmerIdx
        // int k-mer is written into workspace
        void printKmer(int kmerIdx, int kmerSize, char* int2aa);

        // print k amino acids of int k-mer kmer
        void printKmer(const int* kmer, int kmerSize, char* int2aa);
        
        int * powers;

    private:

        int alphabetSize;
        
        int maxKmerSize;

        unsigned int lastKmerIndex;

        unsigned int maxKmerIndex;

        int* workspace;
};
#endif
