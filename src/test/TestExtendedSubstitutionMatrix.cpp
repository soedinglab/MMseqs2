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

const char* binary_name = "test_extendedsubstitutionmatrix";

int main (int argc, const char * argv[])
{
    
    const size_t kmer_size=3;
    
    
    SubstitutionMatrix subMat("/Users/aluucard/Documents/workspace/kClust2/data/blosum30.out",8.0);
    
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("%c\t",subMat.num2aa[i]);
    printf("\n");
//    ReducedMatrix redMat(subMat.probMatrix, subMat.alphabetSize-2);
    
    const int  testSeq[]={1,2,3,1,1,1};
    const int * seq_ptr=&testSeq[0];
    ExtendedSubstitutionMatrix extMat(subMat.subMatrix, kmer_size,subMat.alphabetSize);
    Indexer idx(subMat.alphabetSize,kmer_size);
    
    
    
    std::cout << "Sequence (id 0):\n";
    char* sequence = "AAMICPAEAGRPSLADS";
    std::cout << sequence << "\n\n";
    
    Sequence* s = new Sequence (10000, subMat.aa2num, subMat.num2aa, 0, kmer_size, false);
    s->mapSequence(0,"LALA",sequence);
    
    printf("Normal : ");
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("%c\t",subMat.num2aa[i]);
    printf("\nReduced: ");
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("%c\t",subMat.num2aa[i]);
    printf("\nNormal : ");
    for(int i = 65; i<'Z';i++)
        printf("%d\t",subMat.aa2num[i]);
    printf("\nReduced: ");
    for(int i = 65; i<'Z';i++)
        printf("%d\t",subMat.aa2num[i]);
    
    std::cout << "\nInt reduced sequence:\n";
    for (int i = 0; i < s->L; i++)
        std::cout << s->int_sequence[i] << " ";
    std::cout << "\n";
    
    while(s->hasNextKmer()){
        const int * curr_pos= s->nextKmer();
        printf("kmerpos1: %d\tkmerpos2: %d\n",curr_pos[0],curr_pos[1]);
        unsigned int idx_val=idx.int2index(curr_pos);
        std::cout << "Index:    " <<idx_val << "\n";
        std::cout << "MaxScore: " << extMat.scoreMatrix->score[idx_val*extMat.size]<< "\n";
        
    }
    
    int i = 0;
    
    printf("%d\n",i);
    return 0;
}

