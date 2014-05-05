//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include "Sequence.h"
#include "Indexer.h"
#include "ExtendedSubstitutionMatrix.h"
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "KmerGenerator.h"
#include "BaseMatrix.h"

int main (int argc, const char * argv[])
{

    const size_t kmer_size=5;


    SubstitutionMatrix subMat("../../data/blosum62.out",8.0);
    std::cout << "Subustitution matrix:\n";

    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    //    ReducedMatrix subMat(subMat.probMatrix, 20);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    const int  testSeq[]={1,2,3,1,1,1};
    ExtendedSubstitutionMatrix extMattwo(subMat.subMatrix, 2,subMat.alphabetSize);
    ExtendedSubstitutionMatrix extMatthree(subMat.subMatrix, 3,subMat.alphabetSize);

    Indexer idx(subMat.alphabetSize,kmer_size);
    std::cout << "Sequence (id 0):\n";
    const char* sequence = argv[1];
    std::cout << sequence << "\n\n";
    Sequence* s = new Sequence (10000, subMat.aa2int, subMat.int2aa, 0);
    s->mapSequence(0,"lala",sequence);

    KmerGenerator kmerGen(kmer_size,subMat.alphabetSize,80);

    kmerGen.setDivideStrategy(extMatthree.scoreMatrix, extMattwo.scoreMatrix );
    int* testKmer = new int[kmer_size];
    while(s->hasNextKmer(kmer_size)){
        const int * curr_pos = s->nextKmer(kmer_size);
        printf("kmerpos1: %d\tkmerpos2: %d\n",curr_pos[0],curr_pos[1]);

        unsigned int idx_val=idx.int2index(curr_pos);
        std::cout << "Index:    " <<idx_val << "\n";
//        std::cout << "MaxScore: " << extMattwo.scoreMatrix[idx_val]->back().first<< "\n";
        ScoreMatrix kmer_list= kmerGen.generateKmerList(curr_pos);
        std::cout << "Similar k-mer list size:" << kmer_list.elementSize << "\n\n";

        std::cout << "Similar " << kmer_size << "-mer list for pos 0:\n";
        for (size_t pos = 0; pos < kmer_list.elementSize; pos++){
	    std::cout << "Pos:" << pos << " ";
            std::cout << "Score:" << kmer_list.score[pos]  << " ";
            std::cout << "Index:" << kmer_list.index[pos] << "\n";

            idx.index2int(testKmer, kmer_list.index[pos], kmer_size);
            std::cout << "\t";
            for (size_t i = 0; i < kmer_size; i++)
                std::cout << testKmer[i] << " ";
            std::cout << "\t";
            for (size_t i = 0; i < kmer_size; i++)
                std::cout << subMat.int2aa[testKmer[i]];
            std::cout << "\n";
        }
    }
    return 0;
}

