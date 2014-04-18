#include "BaseMatrix.h"

BaseMatrix::BaseMatrix(){
    this->alphabetSize = 21;
    // init [amino acid <-> int] mappings
    int2aa = new char[alphabetSize];
    int2aa[0] = 'A';
    int2aa[1] = 'R';
    int2aa[2] = 'N';
    int2aa[3] = 'D';
    int2aa[4] = 'C';
    int2aa[5] = 'Q';
    int2aa[6] = 'E';
    int2aa[7] = 'G';
    int2aa[8] = 'H';
    int2aa[9] = 'I';
    int2aa[10] = 'L';
    int2aa[11] = 'K';
    int2aa[12] = 'M';
    int2aa[13] = 'F';
    int2aa[14] = 'P';
    int2aa[15] = 'S';
    int2aa[16] = 'T';
    int2aa[17] = 'W';
    int2aa[18] = 'Y';
    int2aa[19] = 'V';
    int2aa[20] = 'X';

    aa2int = new int['Z'+1];
    for (int i = 0; i <= 'Z'; ++i) aa2int[i]=-1;
    for (int i = 0; i < alphabetSize; ++i){
        aa2int[(int)int2aa[i]] = i;
    }

    // init the background probabilities, joint probability and scoring matrices with zeros
    pBack = new double[alphabetSize];
    probMatrix = new double*[alphabetSize];
    subMatrix = new short*[alphabetSize];

    for (int i = 0; i < alphabetSize; i++){
        pBack[i] = 0.0;
        probMatrix[i] = new double[alphabetSize];
        subMatrix[i] = new short[alphabetSize];
        for (int j = 0; j < alphabetSize; j++){
            probMatrix[i][j] = 0.0;
            subMatrix[i][j] = 0;
        }
    }
}

BaseMatrix::~BaseMatrix(){
    delete[] int2aa;
    delete[] aa2int;
    delete[] pBack;
    for (int i = 0; i < alphabetSize; i++){
        delete[] probMatrix[i];
        delete[] subMatrix[i];
    }
    delete[] probMatrix;
    delete[] subMatrix;
}

void BaseMatrix::print(short** matrix, char* int2aa, int size){
    std::cout << "\n";
    short avg = 0;
    printf("%4c ", ' ');
    for (int i = 0; i < size; i++)
        printf("%4c ", int2aa[i]);
    std::cout << "\n";
    for (int i = 0; i < size; i++){
        printf("%4c ", int2aa[i]);
        for (int j = 0; j < size; j++){
            printf("%4d ", matrix[i][j]);
            avg += matrix[i][j];
        }
        std::cout << "\n";
    }
    std::cout << ((double)avg/(double)(size*size)) << "\n";
}

void BaseMatrix::print(double** matrix, char* int2aa, int size){
    std::cout << "\n";
    double avg = 0.0;
    printf("%7c ", ' ');
    for (int i = 0; i < size; i++)
        printf("%7c ", int2aa[i]);
    std::cout << "\n";
    for (int i = 0; i < size; i++){
        printf("%7c ", int2aa[i]);
        for (int j = 0; j < size; j++){
            printf("%7.4f ", matrix[i][j]);
            avg += matrix[i][j];
        }
        std::cout << "\n";
    }
    std::cout << (avg/(double)(size*size)) << "\n";
}

void BaseMatrix::generateSubMatrix(double ** probMatrix, double ** subMatrix, int size, double bitFactor, double scoringBias){

    // calculate background distribution for the amino acids
    double pBack[size];
    for (int i = 0; i < size; i++){
        pBack[i] = 0;
        for (int j = 0; j < size; j++){
            pBack[i] += probMatrix[i][j];
        }
    }
    // calculate the substitution matrix
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            subMatrix[i][j] = bitFactor * _log2(probMatrix[i][j]/(pBack[i]*pBack[j])) + scoringBias;
        }
    }

    subMatrix[size-1][size-1] = 0.0;
}

void BaseMatrix::generateSubMatrix(double ** probMatrix, short ** subMatrix, int size, double bitFactor, double scoringBias){
    double** sm = new double* [size];
    for (int i = 0; i < size; i++)
        sm[i] = new double[size];

    generateSubMatrix(probMatrix, sm, size, bitFactor, scoringBias);

    // convert to short data type matrix
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            subMatrix[i][j] = (short) floor (sm[i][j] + 0.5);
        }
    }

    for (int i = 0; i < size; i++)
        delete[] sm[i];
    delete[] sm;
}


