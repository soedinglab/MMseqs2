#include <simd/simd.h>
#include <climits> 
#include "BaseMatrix.h"

#include "Debug.h"
#include "Util.h"
#include "MathUtil.h"
#include "Sequence.h"

const double BaseMatrix::ANY_BACK = 1E-5;

BaseMatrix::BaseMatrix(){
    this->alphabetSize = 21;
    // init [amino acid <-> int] mappings

    int2aa = new char[32];
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
    // needed for profile states
    int2aa[21] = 'Z';
    int2aa[22] = '[';
    int2aa[23] = '\\';
    int2aa[24] = ']';
    int2aa[25] = '^';
    int2aa[26] = '_';
    int2aa[27] = '`';
    int2aa[28] = 'a';
    int2aa[29] = 'b';
    int2aa[30] = 'c';
    int2aa[31] = 'd';


    aa2int = new int[UCHAR_MAX];
    for (int i = 0; i < UCHAR_MAX; ++i) aa2int[i]=-1;
    for (int i = 0; i < alphabetSize; ++i){
        aa2int[(int)int2aa[i]] = i;
    }

    // init the background probabilities, joint probability and scoring matrices with zeros
    pBack = new double[alphabetSize];
    probMatrix = new double*[alphabetSize];
    subMatrix = new short*[alphabetSize];
    subMatrix2Bit = new short*[alphabetSize];
    subMatrixPseudoCounts = new float*[alphabetSize];

    for (int i = 0; i < alphabetSize; i++){
        pBack[i] = 0.0;
        probMatrix[i] = new double[alphabetSize];
        subMatrix[i] = new short[alphabetSize];
        subMatrix2Bit[i] = new short[alphabetSize];
        subMatrixPseudoCounts[i] =  (float *) malloc_simd_float(alphabetSize * sizeof(float));
        for (int j = 0; j < alphabetSize; j++){
            probMatrix[i][j] = 0.0;
            subMatrix2Bit[i][j] = 0.0;
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
        delete[] subMatrix2Bit[i];
        free(subMatrixPseudoCounts[i]);
    }
    delete[] probMatrix;
    delete[] subMatrix2Bit;
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

void BaseMatrix::computeBackground(double ** probMat, double * pBack, int alphabetSize, bool containsX){
    for(int i = 0; i < alphabetSize; i++){
        pBack[i] = 0;
        for(int j = 0; j < alphabetSize; j++){
            pBack[i]+=probMat[i][j];
        }
    }
    if(containsX){
        pBack[alphabetSize-1] = ANY_BACK;
    }

}

void BaseMatrix::generateSubMatrix(double ** probMatrix, double ** subMatrix, float ** subMatrixPseudoCounts,
                                   int size, bool containsX) {

    // calculate background distribution for the amino acids
    double *pBack = new double[size];
    computeBackground(probMatrix, pBack, size, containsX);

    //Precompute matrix R for amino acid pseudocounts:
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            subMatrixPseudoCounts[i][j] = probMatrix[i][j] / (pBack[j]); //subMatrixPseudoCounts[a][b]=P(a|b)
        }
    }

    // calculate the substitution matrix
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            subMatrix[i][j] = std::log2(probMatrix[i][j] / (pBack[i] * pBack[j]));
            //printf("%3.4f %3.4f %3.4f \n",  probMatrix[i][j], pBack[i], pBack[j] );

        }
    }
    delete[] pBack;
//    subMatrix[size - 1][size - 1] = -.7;
//    for (int i = 0; i < size; i++) {
//        subMatrix[size - 1][i] = -.7;
//        subMatrix[i][size - 1] = -.7;
//    }
}
void BaseMatrix::generateSubMatrix(double ** probMatrix, float ** subMatrixPseudoCounts, short ** subMatrix,
                                   short **subMatrix2Bit, int size, bool containsX, double bitFactor, double scoringBias){
    double** sm = new double* [size];
    for (int i = 0; i < size; i++)
        sm[i] = new double[size];

    generateSubMatrix(probMatrix, sm, subMatrixPseudoCounts, size, containsX);

    // convert to short data type matrix
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            double pValNBitScale = (bitFactor * sm[i][j] + scoringBias);
            subMatrix[i][j] = (pValNBitScale < 0.0) ? pValNBitScale - 0.5 : pValNBitScale + 0.5;
            double pVal2BitScale = (sm[i][j] * 2.0 + scoringBias);
            subMatrix2Bit[i][j] = (pVal2BitScale < 0.0) ? pVal2BitScale - 0.5 : pVal2BitScale + 0.5;
        }
    }

    for (int i = 0; i < size; i++)
        delete[] sm[i];
    delete[] sm;
}

std::string BaseMatrix::getMatrixName() {
    return matrixName;
}

double BaseMatrix::getBackgroundProb(size_t aa_index)  {
    Debug(Debug::ERROR) << "getBackground is not Impl. for this type of Matrix \n";
    EXIT(EXIT_FAILURE);
}

