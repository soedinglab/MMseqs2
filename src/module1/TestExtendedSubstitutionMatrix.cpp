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

int main (int argc, const char * argv[])
{
    
    const size_t kmer_size=3;
    
    
    SubstitutionMatrix subMat("/cluster/user/maria/kClust2/data/blosum30.out",20);
    
    for(int i = 0; i<subMat.ALPHABET_SIZE;i++)
        printf("%c\t",subMat.int2aa[i]);
    printf("\n");
    ReducedMatrix redMat(subMat.probMatrix,
                        subMat.aa2int,subMat.int2aa,subMat.ALPHABET_SIZE,subMat.ALPHABET_SIZE-2);
    
    const int  testSeq[]={1,2,3,1,1,1};
    const int * seq_ptr=&testSeq[0];
    ExtendedSubstitutionMatrix extMat(redMat.reduced_Matrix, kmer_size,redMat.reduced_alphabet_size);
    Indexer idx(redMat.reduced_alphabet_size,kmer_size);
    
    
    
    std::cout << "Sequence (id 0):\n";
    char* sequence = "AAMICPAEAGRPSLADS";
    std::cout << sequence << "\n\n";
    
    Sequence* s = new Sequence (10000, redMat.reduced_aa2int, redMat.reduced_int2aa);
    s->setId(0);
    s->mapSequence(sequence);
    
    printf("Normal : ");
    for(int i = 0; i<subMat.ALPHABET_SIZE;i++)
        printf("%c\t",subMat.int2aa[i]);
    printf("\nReduced: ");
    for(int i = 0; i<subMat.ALPHABET_SIZE;i++)
        printf("%c\t",redMat.reduced_int2aa[i]);
    printf("\nNormal : ");
    for(int i = 65; i<'Z';i++)
        printf("%d\t",subMat.aa2int[i]); 
    printf("\nReduced: ");
    for(int i = 65; i<'Z';i++)
        printf("%d\t",redMat.reduced_aa2int[i]); 
    
    std::cout << "\nInt reduced sequence:\n";
    for (int i = 0; i < s->L; i++)
        std::cout << s->int_sequence[i] << " ";
    std::cout << "\n";
    
    while(s->hasNextKmer(kmer_size)){
        const int * curr_pos= s->nextKmer(kmer_size);
        printf("kmerpos1: %d\tkmerpos2: %d\n",curr_pos[0],curr_pos[1]);
        size_t idx_val=idx.int2index(curr_pos);
        std::cout << "Index:    " <<idx_val << "\n";
        std::cout << "MaxScore: " << extMat.scoreMatrix[idx_val]->at(0).first<< "\n";
        
    }
    
    int i = 0;
    
    printf("%d\n",i);
    return 0;
}

