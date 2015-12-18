#include <string>
#include "BaseMatrix.h"

BaseMatrix::BaseMatrix(){
    this->alphabetSize = 21;
    // init [amino acid <-> int] mappings
    int2aa = new char[alphabetSize];
    // A C D E F G	H I	K L M N P Q R S T V W Y
    int2aa[0] = 'A';
    int2aa[1] = 'C';
    int2aa[2] = 'D';
    int2aa[3] = 'E';
    int2aa[4] = 'F';
    int2aa[5] = 'G';
    int2aa[6] = 'H';
    int2aa[7] = 'I';
    int2aa[8] = 'K';
    int2aa[9] = 'L';
    int2aa[10] = 'M';
    int2aa[11] = 'N';
    int2aa[12] = 'P';
    int2aa[13] = 'Q';
    int2aa[14] = 'R';
    int2aa[15] = 'S';
    int2aa[16] = 'T';
    int2aa[17] = 'V';
    int2aa[18] = 'W';
    int2aa[19] = 'Y';
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
    subMatrixDouble = new double*[alphabetSize];
    subMatrixPseudoCounts = new float*[alphabetSize];

    for (int i = 0; i < alphabetSize; i++){
        pBack[i] = 0.0;
        probMatrix[i] = new double[alphabetSize];
        subMatrix[i] = new short[alphabetSize];
        subMatrixDouble[i] = new double[alphabetSize];
        subMatrixPseudoCounts[i] = new float[alphabetSize];
        for (int j = 0; j < alphabetSize; j++){
            probMatrix[i][j] = 0.0;
            subMatrixDouble[i][j] = 0.0;
            subMatrix[i][j] = 0;
            subMatrixPseudoCounts[i][j] = 0.0;
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
        delete[] subMatrixPseudoCounts[i];
    }
    delete[] probMatrix;
    delete[] subMatrixPseudoCounts;
    delete[] subMatrix;
}

void BaseMatrix::print(short** matrix, char* int2aa, int size){
    std::cout << "\n";
    short avg = 0;
    printf("     ");
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

void BaseMatrix::generateSubMatrix(double ** probMatrix, double ** subMatrix, float ** subMatrixPseudoCounts,
                                   int size, double bitFactor, double scoringBias){

    // calculate background distribution for the amino acids
    double pBack[size];
    for (int i = 0; i < size; i++){
        pBack[i] = 0;
        for (int j = 0; j < size; j++){
            pBack[i] += probMatrix[i][j];
        }
    }
    //Precompute matrix R for amino acid pseudocounts:
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            subMatrixPseudoCounts[i][j] = probMatrix[i][j]/(pBack[j]); //subMatrixPseudoCounts[a][b]=P(a|b)
        }
    }

    // calculate the substitution matrix
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            subMatrix[i][j] = bitFactor * _log2(probMatrix[i][j]/(pBack[i]*pBack[j])) + scoringBias;

        }
    }

    //subMatrix[size-1][size-1] = 0.0;
}

void BaseMatrix::generateSubMatrix(double ** probMatrix, float ** subMatrixPseudoCounts, short ** subMatrix,
                                   double ** subMatrixDouble, int size, double bitFactor, double scoringBias){
    double** sm = new double* [size];
    for (int i = 0; i < size; i++)
        sm[i] = new double[size];

    generateSubMatrix(probMatrix, sm, subMatrixPseudoCounts, size, bitFactor, scoringBias);

    // convert to short data type matrix
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            subMatrix[i][j] = (short) floor (sm[i][j] + 0.5);
            subMatrixDouble[i][j] = sm[i][j];
        }
    }

    for (int i = 0; i < size; i++)
        delete[] sm[i];
    delete[] sm;
}

std::string BaseMatrix::getMatrixName() {
    return matrixName;
}
