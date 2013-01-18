//
//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#include "ReducedMatrix.h"
#include "SubstitutionMatrix.h"

int main (int argc, const char * argv[])
{
    
    SubstitutionMatrix subMat("../../data/blosum30.out",20);
    
    for(int i = 0; i<subMat.ALPHABET_SIZE;i++)
        printf("(%d)%c\t",i,subMat.int2aa[i]);
    printf("\n");
    ReducedMatrix redMat(subMat.probMatrix,
                        subMat.aa2int,subMat.int2aa,subMat.ALPHABET_SIZE,16);
    std::cout << "\n";
    printf("Normal alphabet : ");
    for(int i = 0; i<subMat.ALPHABET_SIZE;i++)
        printf("%c\t",subMat.int2aa[i]);
    printf("\nReduced alphabet: ");
    for(int i = 0; i<subMat.ALPHABET_SIZE;i++)
        printf("%c\t",redMat.reduced_int2aa[i]);
    std::cout << "\nReduced alphabet size: " << redMat.reduced_alphabet_size << "\n";

    std::cout << "aa2int: \n";
    for (char c = 'A'; c <= 'Z'; c++)
        printf("%c%3d\t", c, subMat.aa2int[c]);
    std::cout << "\n";

    std::cout << "reduced aa2int:\n";
    for (char c = 'A'; c <= 'Z'; c++)
        printf("%c%3d\t", c, redMat.reduced_aa2int[c]);
    std::cout << "\n";

    printf("\n\nOriginal substitution matrix:\n");
    printf("   \t");
    for (int i = 0; i < subMat.ALPHABET_SIZE; i++)
        printf("%3c\t",subMat.int2aa[i]);
    std::cout << "\n";
    for (int i = 0; i < subMat.ALPHABET_SIZE; i++){
        printf("%3c\t",subMat.int2aa[i]);
        for (int j = 0; j < subMat.ALPHABET_SIZE; j++){
            printf("%3d\t", subMat.scMatrix[i][j]);
        }
        std::cout << "\n";
    }

    printf("\n\nReduced substitution matrix:\n");
    printf("   \t");
    for(size_t i = 0; i<redMat.reduced_alphabet_size;i++)
        printf("%3c\t",redMat.reduced_alphabet->at(i));
    std::cout << "\n";
    for (int i = 0; i < redMat.reduced_alphabet_size; i++){
        printf("%3c\t",redMat.reduced_alphabet->at(i));
        for (int j = 0; j < redMat.reduced_alphabet_size; j++){
            printf("%3d\t", redMat.reduced_Matrix[i][j]);
        }
        std::cout << "\n";
    }

    return 0;
}

