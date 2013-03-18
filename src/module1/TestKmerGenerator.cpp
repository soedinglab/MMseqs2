//
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
    
    
    SubstitutionMatrix subMat("/Users/aluucard/Documents/workspace/kClust2/data/blosum62.out");
    std::cout << "Subustitution matrix:\n";
 //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    std::cout << "\n";

    std::cout << "subMatrix:\n";
    ReducedMatrix redMat(subMat.probMatrix, 20);
 //   BaseMatrix::print(redMat.subMatrix, redMat.alphabetSize);
    std::cout << "\n";
    
    const int  testSeq[]={1,2,3,1,1,1};
    ExtendedSubstitutionMatrix extMattwo(redMat.subMatrix, 2,redMat.alphabetSize);
    ExtendedSubstitutionMatrix extMatthree(redMat.subMatrix, 3,redMat.alphabetSize);

    Indexer idx(redMat.alphabetSize,kmer_size);
    
    
    
    std::cout << "Sequence (id 0):\n";
    char* sequence = "AAVIDE";
    std::cout << sequence << "\n\n";
    
    Sequence* s = new Sequence (10000, redMat.aa2int, redMat.int2aa);
    s->mapSequence(sequence);
    
    printf("Normal alphabet : ");
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("%c\t",subMat.int2aa[i]);
    
    printf("\nReduced alphabet: ");
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("%c\t",redMat.int2aa[i]);
    
    printf("\nNormal int code: ");
    for(int i = 'A'; i<'Z';i++)
        printf("%d\t",subMat.aa2int[i]); 
    
    printf("\nReduced int code: ");
    for(int i = 'A'; i<'Z';i++)
        printf("%d\t",redMat.aa2int[i]); 
    
    std::cout << "\nInt reduced sequence:\n";
    for (int i = 0; i < s->L; i++)
        std::cout << s->int_sequence[i] << " ";
    std::cout << "\nChar reduced sequence:\n";
    for (int i = 0; i < s->L; i++)
        std::cout << redMat.int2aa[s->int_sequence[i]] << " ";
    std::cout << "\n";
    
    KmerGenerator kmerGen(kmer_size,redMat.alphabetSize,10, 
                          &extMatthree,&extMattwo );
    
    int* testKmer = new int[kmer_size];
    while(s->hasNextKmer(kmer_size)){
        const int * curr_pos = s->nextKmer(kmer_size);
        printf("kmerpos1: %d\tkmerpos2: %d\n",curr_pos[0],curr_pos[1]);
        
        unsigned int idx_val=idx.int2index(curr_pos);
        std::cout << "Index:    " <<idx_val << "\n";
//        std::cout << "MaxScore: " << extMattwo.scoreMatrix[idx_val]->back().first<< "\n";
        
        KmerGeneratorResult kmer_list= kmerGen.generateKmerList(curr_pos);
        
        std::pair<short,unsigned int> ** retList=kmer_list.scoreKmerList;
        
        std::cout << "Similar k-mer list size:" << kmer_list.count << "\n\n";

        std::cout << "Similar " << kmer_size << "-mer list for pos 0:\n";
        for (int pos = 0; pos < kmer_list.count; pos++){
            std::pair<short,unsigned int> * result = retList[pos];
            std::cout << "Score:" << result->first << "\n";
            std::cout << "Index:" << result->second << "\n";

            idx.index2int(testKmer, result->second, kmer_size);
            std::cout << "\t";
            for (int i = 0; i < kmer_size; i++)
                std::cout << testKmer[i] << " ";
            std::cout << "\t";
            for (int i = 0; i < kmer_size; i++)
                std::cout << redMat.int2aa[testKmer[i]];
            std::cout << "\n";
        }
    }
    
    
    return 0;
}

