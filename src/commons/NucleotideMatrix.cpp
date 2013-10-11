#include "NucleotideMatrix.h"

NucleotideMatrix::NucleotideMatrix(){

    this->alphabetSize = 5;
    delete[] int2aa;
    delete[] aa2int;

    int2aa = new char[alphabetSize];
    int2aa[0] = 'a';
    int2aa[1] = 'c';
    int2aa[2] = 'g';
    int2aa[3] = 't';
    int2aa[4] = 'n';

    aa2int = new int['z'+1];
    for (int i = 0; i <= 'x'; ++i) aa2int[i]=-1;
    for (int i = 0; i < alphabetSize; ++i){
        aa2int[(int)int2aa[i]] = i;
    }

    subMatrix = new short*[alphabetSize];
    for (int i = 0; i < alphabetSize; i++)
        subMatrix[i] = new short[alphabetSize];

    for (int i = 0; i < alphabetSize-1; i++){
        for (int j = 0; j < alphabetSize-1; j++){
            if (i == j)
                subMatrix[i][j] = 3;
            else
                subMatrix[i][j] = -2;
        }
    }

    for (int i = 0; i < alphabetSize; i++){
        subMatrix[alphabetSize-1][i] = 0;
        subMatrix[i][alphabetSize-1] = 0;
    }

}

NucleotideMatrix::~NucleotideMatrix(){
}

