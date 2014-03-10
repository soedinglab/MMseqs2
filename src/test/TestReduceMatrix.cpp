//
//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#include "BaseMatrix.h"
#include "ReducedMatrix.h"
#include "SubstitutionMatrix.h"

int main (int argc, const char * argv[])
{
    
    SubstitutionMatrix subMat("../../data/blosum62.out",8.0);
    
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("(%d)%c\t",i,subMat.int2aa[i]);
    printf("\n");
    ReducedMatrix redMat(subMat.probMatrix, 5);
    std::cout << "\n";
    printf("Normal alphabet : ");
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("%c\t",subMat.int2aa[i]);
    printf("\nReduced alphabet: ");
    for(int i = 0; i<redMat.alphabetSize;i++)
        printf("%c\t",redMat.int2aa[i]);
    std::cout << "\nReduced alphabet size: " << redMat.alphabetSize << "\n";

    std::cout << "aa2int: \n";
    for (char c = 'A'; c <= 'Z'; c++)
        printf("%c%3d\t", c, subMat.aa2int[c]);
    std::cout << "\n";

    std::cout << "int2aa: \n";
    for (int i = 0; i < subMat.alphabetSize; i++)
        printf("%d%3c\t", i, subMat.int2aa[i]);
    std::cout << "\n";

    std::cout << "reduced aa2int:\n";
    for (char c = 'A'; c <= 'Z'; c++)
        printf("%c%3d\t", c, redMat.aa2int[c]);
    std::cout << "\n";

    std::cout << "reduced int2aa: \n";
    for (int i = 0; i < redMat.alphabetSize; i++)
        printf("%d%3c\t", i, redMat.int2aa[i]);
    std::cout << "\n";

    printf("\n\nOriginal substitution matrix:\n");
//    subMat.print();

    printf("\n\nReduced substitution matrix:\n");
//    redMat.print();

    return 0;
}

