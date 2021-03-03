#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include <string>

class BaseMatrix{
public:
    BaseMatrix();

    virtual ~BaseMatrix();

    /*contains amino acid to int mapping*/
    unsigned char* aa2num;

    /*contains int to amino acid mapping*/
    char* num2aa;

    /* size of alphabet*/
    int alphabetSize;

    // allocated alphabet size
    int allocatedAlphabetSize;

    /* matrixData */
    std::string matrixData;

    // substitution matrix
    short** subMatrix;
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

    // score bias for simpleGotoh
    float scoreBias;

    // print the substitution matrix
    static void print(short** matrix, char* num2aa, int size);

    static void print(double** matrix, char* num2aa, int size);

    void initMatrixMemory(int alphSize);

    // generate the substitution matrix given the probability matrix, background probabilities and the alphabet size
    static void generateSubMatrix(double ** probMatrix, double ** subMatrix, float ** subMatrixPseudoCounts, int size, bool containsX);

    // generate a short data type substitution matrix
    void generateSubMatrix(double ** probMatrix, float ** subMatrixPseudoCounts, short ** subMatrix, int size, bool containsX, double bitFactor = 1.0, double scoringBias = 0.0);

    virtual double getBackgroundProb(size_t aa_index);

    virtual void setupLetterMapping() {};

    virtual float getBitFactor() {return 1.0; }

    std::string getMatrixName();

    std::string matrixName;

    inline double getLambda() {
        return lambda;
    }

    static void computeBackground(double **probMat, double *pBack, int alphabetSize, bool containsX);

    static size_t memorySize(std::string & matrixName , std::string & matrixData);
    static std::pair<std::string, std::string> unserialize(const char * data);
    static char * serialize(std::string &matrixName, std::string &matrixData );
    static std::string unserializeName(const char * data);
};


class ProbabilityMatrix {
public:
    ProbabilityMatrix(BaseMatrix &matrix) : alphabetSize(matrix.alphabetSize) {
        probMatrix = new double*[matrix.alphabetSize];
        probMatrixPointers = new const double*[matrix.alphabetSize];
        std::fill_n(hardMaskTable, 256, matrix.aa2num[static_cast<int>('X')]);
        for (int i = 0; i < matrix.alphabetSize; ++i) {
            probMatrix[i] = new double[matrix.alphabetSize];
            probMatrixPointers[i] = probMatrix[i];
            for (int j = 0; j < matrix.alphabetSize; ++j) {
                probMatrix[i][j] = matrix.probMatrix[i][j] / (matrix.pBack[i] * matrix.pBack[j]);
            }
        }
    }
    ~ProbabilityMatrix() {
        for (int i = 0; i < alphabetSize; ++i) {
            delete[] probMatrix[i];
        }
        delete[] probMatrix;
        delete[] probMatrixPointers;
    }

    char hardMaskTable[256];
    const double **probMatrixPointers;

private:
    const int alphabetSize;
    double **probMatrix;

};
#endif
