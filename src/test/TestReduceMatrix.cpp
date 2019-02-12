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
#include "Parameters.h"

const char* binary_name = "test_reducematrix";

int main (int, const char**) {
    const int reductionAlphabetSize = 17;
    Parameters& par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.c_str(), 2.0, 0);
    subMat.print(subMat.subMatrix, subMat.int2aa,21);
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("(%d, %c) ",i,subMat.int2aa[i]);
    printf("\n");
    ReducedMatrix redMat(subMat.probMatrix, subMat.subMatrixPseudoCounts, subMat.aa2int, subMat.int2aa, subMat.alphabetSize, reductionAlphabetSize, subMat.getBitFactor());
    std::cout << "\n";
    printf("Normal alphabet : ");
    for(int i = 0; i<subMat.alphabetSize;i++)
        printf("(%c) ",subMat.int2aa[i]);
    printf("\nReduced alphabet: ");
    for(int i = 0; i<redMat.alphabetSize;i++)
        printf("(%c) ",redMat.int2aa[i]);
    std::cout << "\nReduced alphabet size: " << redMat.alphabetSize << "\n";

    std::cout << "aa2int: \n";
    for (char c = 'A'; c <= 'Z'; c++)
        printf("%c%3d\t", c, subMat.aa2int[(int)c]);
    std::cout << "\n";

    std::cout << "int2aa: \n";
    for (int i = 0; i < subMat.alphabetSize; i++)
        printf("%d%3c\t", i, subMat.int2aa[i]);
    std::cout << "\n";

    std::cout << "reduced aa2int:\n";
    for (char c = 'A'; c <= 'Z'; c++)
        printf("%c%3d\t", c, redMat.aa2int[(int)c]);
    std::cout << "\n";

    std::cout << "reduced int2aa: \n";
    for (int i = 0; i < redMat.alphabetSize; i++)
        printf("%d%3c\t", i, redMat.int2aa[i]);
    std::cout << "\n";

    printf("\n\nOriginal substitution matrix:\n");
    subMat.print(subMat.subMatrix, subMat.int2aa,21);
    subMat.print(subMat.probMatrix, subMat.int2aa,21);

    printf("\n\nReduced substitution matrix:\n");
    subMat.print(redMat.subMatrix, redMat.int2aa,reductionAlphabetSize);
    subMat.print(redMat.probMatrix, redMat.int2aa,reductionAlphabetSize);


    return 0;
}

