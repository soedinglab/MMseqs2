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
#include "Parameters.h"

const char* binary_name = "test_kmergenerator";

int main (int, const char**) {
    const size_t kmer_size=6;

    Parameters& par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.values.aminoacid().c_str(), 8.0, 0);
    std::cout << "Subustitution matrix:\n";

    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    //    ReducedMatrix subMat(subMat.probMatrix, 20);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout <<  std::endl;
    std::cout << "ExtSupMatrix:"<< std::endl;

    ScoreMatrix extMattwo = ExtendedSubstitutionMatrix::calcScoreMatrix(subMat, 2);
    ScoreMatrix extMatthree = ExtendedSubstitutionMatrix::calcScoreMatrix(subMat, 3);

    Indexer idx(subMat.alphabetSize,kmer_size);
    std::cout << "Sequence (id 0):\n";
    const char* sequence = "PATWPCLVALG";
    std::cout << sequence << "\n\n";
    Sequence* s = new Sequence(10000, Parameters::DBTYPE_AMINO_ACIDS, &subMat, kmer_size, false, false);
    s->mapSequence(0,0,sequence, strlen(sequence));

    KmerGenerator kmerGen(kmer_size,subMat.alphabetSize,161);
    kmerGen.setDivideStrategy(&extMatthree, &extMattwo);

    size_t * testKmer = new size_t[kmer_size];
    int i = 0; 
    while(s->hasNextKmer()){
        const unsigned char * curr_pos = s->nextKmer();
        printf("Pos1: %d\n", i++);

        unsigned int idx_val=idx.int2index(curr_pos);
        std::cout << "Index:    " <<idx_val << "  ";
        idx.printKmer(idx_val, kmer_size, subMat.num2aa);
        std::cout << std::endl;
//        std::cout << "MaxScore: " << extMattwo.scoreMatrix[idx_val]->back().first<< "\n";
        std::pair<size_t *, size_t > kmer_list= kmerGen.generateKmerList(curr_pos);
        std::cout << "Similar k-mer list size:" << kmer_list.second << "\n\n";

        std::cout << "Similar " << kmer_size << "-mer list for pos 0:\n";
        for (size_t pos = 0; pos < kmer_list.second; pos++){
	    std::cout << "Pos:" << pos << " ";
//            std::cout << "Score:" << kmer_list.score[pos]  << " ";
            std::cout << "Index:" << kmer_list.first[pos] << "\n";

            idx.index2int(testKmer, kmer_list.first[pos], kmer_size);
            std::cout << "\t";
            for (size_t i = 0; i < kmer_size; i++)
                std::cout << testKmer[i] << " ";
            std::cout << "\t";
            for (size_t i = 0; i < kmer_size; i++)
                std::cout << subMat.num2aa[testKmer[i]];
            std::cout << "\n";
        }
    }

    ExtendedSubstitutionMatrix::freeScoreMatrix(extMatthree);
    ExtendedSubstitutionMatrix::freeScoreMatrix(extMattwo);

    return EXIT_SUCCESS;
}

