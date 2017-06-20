#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include <string>

class BaseMatrix{
public:
    BaseMatrix();

    virtual ~BaseMatrix();

    /*contains amino acid to int mapping*/
    int*  aa2int;

    /*contains int to amino acid mapping*/
    char* int2aa;

    /* size of alphabet*/
    int alphabetSize;

    // substitution matrix
    short** subMatrix;

    //  substitution matrix in double
    short **subMatrix2Bit;

    // substitution matrix for pseudocounts
    float** subMatrixPseudoCounts;

    // joint probability matrix
    double** probMatrix;

    // scaling factor computed by cacode
    double lambda;

    // background probabilities of the amino acids
    double* pBack;

    // background for any state
    static const double ANY_BACK;

    // print the substitution matrix
    static void print(short** matrix, char* int2aa, int size);

    static void print(double** matrix, char* int2aa, int size);

    // generate the substitution matrix given the probability matrix, background probabilities and the alphabet size
    static void generateSubMatrix(double ** probMatrix, double ** subMatrix, float ** subMatrixPseudoCounts, int size, bool containsX);

    // generate a short data type substitution matrix
    static void generateSubMatrix(double ** probMatrix, float ** subMatrixPseudoCounts, short ** subMatrix, short **subMatrix2Bit, int size, bool containsX, double bitFactor = 1.0, double scoringBias = 0.0);

    virtual double getBackgroundProb(size_t aa_index);

    virtual void setupLetterMapping() {};

    virtual float getBitFactor() {return 1.0; }

    std::string getMatrixName();

    std::string matrixName;

    inline double getLambda() {
        return lambda;
    }

    static void computeBackground(double **probMat, double *pBack, int alphabetSize, bool containsX);
};
#endif
