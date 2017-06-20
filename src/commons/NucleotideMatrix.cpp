#include "NucleotideMatrix.h"

NucleotideMatrix::NucleotideMatrix(){
    this->alphabetSize = 5;

    setupLetterMapping();
    matrixName = "NUCL";


    for (int i = 0; i < alphabetSize-1; i++){
        for (int j = 0; j < alphabetSize-1; j++){
            if (i == j)
                subMatrix[i][j] = 3;
            else
                subMatrix[i][j] = -2;
        }
    }

    for (int i = 0; i < alphabetSize; i++){
        subMatrix[alphabetSize-1][i] = -1;
        subMatrix[i][alphabetSize-1] = -1;
    }
}


void NucleotideMatrix::setupLetterMapping(){

}


NucleotideMatrix::~NucleotideMatrix(){
}

