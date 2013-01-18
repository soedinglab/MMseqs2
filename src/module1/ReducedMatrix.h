#ifndef ReducedMatrix_H
#define ReducedMatrix_H
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

class ReducedMatrix{
    public: 
        ReducedMatrix(double **probMatrix,int * aa2int,char * int2aa,
                size_t old_alphabet_size, size_t reduced_alphabet_size);
        ~ReducedMatrix();
        /*contains reduced matrix*/
        short ** reduced_Matrix;
        /*mapping aa (orig. alphabet) -> int code of the representative amino acid*/
        int*   reduced_aa2int;
        /*mapping int code (orig. alphabet) -> the representative amino acid char*/
        char*  reduced_int2aa;
        /* size of reduced alphabet*/
        size_t reduced_alphabet_size;

        std::vector<char>* reduced_alphabet;

    private: 
        /*contains original amino acid to int mapping*/
        int*   aa2int;
        /*contains original int to amino acid mapping*/
        char*  int2aa;
        /* size of original alphabet*/
        size_t alphabet_size;

        inline double _log2 (double x) { return log10(x)/0.301029996; }   
        /* The function adds two rows of a given m x n input matrix and produces a
         * m-1 x n matrix. row1 and row2 (where row1 < row2) are the rows that are
         * to be added. row1 is replaced by row1 + row2 in the output matrix.
         */
        void addTwoRows(double ** input, double ** output, size_t numRows , size_t numCols, size_t row1, size_t row2 );
        /* The function adds two columns of a given m x n input matrix and produces a
         * m x n-1 matrix. col1 and col2 (where col1 < col2) are the columns that are
         * to be added. col1 is replaced by col1 + col2 in the output matrix
         */
        void addTwoColumns(double ** input, double ** output, size_t numRows , size_t numCols, size_t col1, size_t col2 );
        /* This function calculates the independent probability for all bases so that we don't have
         * to calculate them every time while we are generating the substitution matrix.
         */
        void genProbBaseArray(double ** pmatrix, double * prob, size_t numRows, size_t numCols);

        /* Copy from array to array */
        void copyMatrix(double ** input,double ** output, size_t numRows, size_t numCols);
        /* This function generates the corresponding substitution matrix given the corresponding probability
         * matrix.
         */
        void generateSubMatrix(double ** pmatrix, double ** subMatrix, size_t numRows, size_t numCols);
        /* This function generates the corresponding substitution matrix given the corresponding probability
         * matrix.
         */
        void generateBiasedSubMatrix(double ** pmatrix, short ** subMatrix,
                double bitFactor, double scoringBias,
                size_t numRows, size_t numCols);
        /* The function coupleBases takes in an input array with numRows and numCols
         * containing meaningful information. This couples the bases base1 and base2
         * by adding the rows and columns corresponding to it (given base1 < base2).
         */
        void coupleBases(double ** input, double ** output, size_t numRows, size_t numCols, size_t base1, size_t base2);
        /* This function calculates the mutual information of a given substitution
         * matrix (subMatrix) given that its Probability matrix (pMatrix) is also known.
         * numRows and numCols are the number of rows and columns respectively in the
         * matrices that contain meaningful information.
         *
         * Mutual information measures the information that two random variables X and Y share: 
         * it measures how much knowing one of these variables reduces uncertainty about the other.
         * (Knowing the amino acid X - how much we know about the aligned amino acid?).
         */
        double calculateMutualInformation(double ** pMatrix, double ** subMatrix, size_t numRows, size_t numCols);
        /* This function finds the two best bases to couple such that we loose the minimum amount of information.
         * Returns the amount of mutual information in the best pairing.
         */
        std::pair<size_t,size_t> coupleWithBestInfo(double ** pinput, double ** pMatrix, size_t numRows, size_t numCols);
};


#endif  
