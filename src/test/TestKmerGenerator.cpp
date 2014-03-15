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

    const size_t kmer_size=4;

    SubstitutionMatrix subMat("../../data/blosum62.out",8.0);
    std::cout << "Subustitution matrix:\n";
    std::cout << "lala matrix:\n";

 //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
 //   ReducedMatrix redMat(subMat.probMatrix, 20);
 //   BaseMatrix::print(redMat.subMatrix, redMat.alphabetSize);
    std::cout << "\n";
    
    const int  testSeq[]={1,2,3,1,1,1};
    ExtendedSubstitutionMatrix extMattwo(subMat.subMatrix, 2,subMat.alphabetSize);
    ExtendedSubstitutionMatrix extMatthree(subMat.subMatrix,3 ,subMat.alphabetSize);
    

    Indexer idx(subMat.alphabetSize,kmer_size);
    
    
    
    std::cout << "Sequence (id 0):\n";
    char* sequence = (char *) argv[1];
    std::cout << sequence << "\n\n";
    
    Sequence* s = new Sequence (10000, subMat.aa2int, subMat.int2aa,0);
    s->mapSequence(0,"lala",sequence);
    
    printf("Normal alphabet : ");
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("%c\t",subMat.int2aa[i]);
    
    
    printf("\nNormal int code: ");
    for(int i = 'A'; i<'Z';i++)
        printf("%d\t",subMat.aa2int[i]); 
    
    
    KmerGenerator kmerGen(kmer_size,subMat.alphabetSize,10,
                          &extMatthree,&extMattwo );

    int* testKmer = new int[kmer_size];
    while(s->hasNextKmer(kmer_size)){
        const int * curr_pos = s->nextKmer(kmer_size);
        printf("kmerpos1: %d\tkmerpos2: %d\n",curr_pos[0],curr_pos[1]);
        
        unsigned int idx_val=idx.int2index(curr_pos);
        std::cout << "Index:    " <<idx_val << "\n";
//        std::cout << "MaxScore: " << extMattwo.scoreMatrix[idx_val]->back().first<< "\n";
    
        KmerGeneratorResult kmer_list= kmerGen.generateKmerList(curr_pos);
        
        std::cout << "Similar k-mer list size:" << kmer_list.count << "\n\n";

        std::cout << "Similar " << kmer_size << "-mer list for pos 0:\n";
        for (int pos = 0; pos < kmer_list.count; pos++){
            std::cout << "Score:" << kmer_list.score[pos] << "\n";
            std::cout << "Index:" << kmer_list.index[pos] << "\n";
//            idx.printKmer(testKmer, kmer_size, subMat.int2aa);
            idx.index2int(testKmer, kmer_list.index[pos] , kmer_size);
            std::cout << "\t";
            for (int i = 0; i < kmer_size; i++)
                std::cout << testKmer[i] << " ";
            std::cout << "\t";
            for (int i = 0; i < kmer_size; i++)
                std::cout << subMat.int2aa[testKmer[i]];
            std::cout << "\n";
        }
    }
    
    
    return 0;
}

