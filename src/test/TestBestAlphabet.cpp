//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//
#include <iostream>
#include <cmath>
#include "SubstitutionMatrix.h"
#include "ReducedMatrix.h"
#include "BaseMatrix.h"
#include "Parameters.h"
const char* binary_name = "test_bestalphabet";

int main (int, const char**) {
    SubstitutionMatrix subMat("blosum62.out", 2.0, -0.0f);
    std::cout << "Subustitution matrix:\n";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.int2aa,subMat.alphabetSize);
    for(size_t aa1 = 0; aa1 < subMat.alphabetSize; aa1++) {
        for (size_t aa2 = 0; aa2 < subMat.alphabetSize; aa2++) {
            printf("%.4f\t", subMat.probMatrix[aa1][aa2]);
        }
        printf("\n");
    }

    for(size_t aa1 = 0; aa1 < subMat.alphabetSize; aa1++) {
        printf("%.5f\t", subMat.pBack[aa1]);
    }
    printf("\n");

    return 0;
}